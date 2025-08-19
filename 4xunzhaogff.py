#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO  # 用于解析 FASTA 文件

# ========== 配置路径 ==========
input_gff    = "data/ninanjie/ninanjie1/ninanjie1.gff"
input_fna    = "data/ninanjie/ninanjie1/ninanjie1.fna"
input_base   = Path(input_gff).parent.parent
source_dir   = input_base / "output"
output_base  = input_base.parent / "ninanjieworkflow3"
splitgff_dir = output_base / "splitgff"
gffindex_dir = output_base / "gffindex"
annotated_dir = output_base / "sp18"  # 新增注释输出目录
fa_pattern   = str(input_base) + "workflow2/*.fa"

# ========== 确保目录存在 ==========
output_base.mkdir(parents=True, exist_ok=True)
splitgff_dir.mkdir(parents=True, exist_ok=True)
gffindex_dir.mkdir(parents=True, exist_ok=True)
annotated_dir.mkdir(parents=True, exist_ok=True)

# ========== 工具函数 ==========
def format_chr_name(text):
    return re.sub(r'chr(\d)(?=\D|$)', r'chr0\1', text)

def process_file(filepath):
    content = Path(filepath).read_text()
    Path(filepath).write_text(format_chr_name(content))

def unify_chr_numbers():
    print("🔍 统一 chr 编号格式...")
    for bed_dir in glob.glob(str(output_base / "*.bed")):
        fa = Path(bed_dir) / "yuangene.fa"
        if fa.exists():
            process_file(str(fa))
        gff3 = Path(bed_dir) / "yuangene.gff3"
        if gff3.exists():
            process_file(str(gff3))
    print("✅ chr 编号格式统一完成")

def extract_prefix(filename):
    m = re.match(r'^([A-Za-z0-9]+)', filename)
    return m.group(1) if m else "split"

def split_and_rename_gff(file_path, output_dir, prefix=None):
    if prefix is None:
        prefix = extract_prefix(Path(file_path).name)
    lines = Path(file_path).read_text().splitlines()
    groups = defaultdict(list)
    for L in lines:
        if L.startswith('#'):
            continue
        chrom = L.split('\t',1)[0]
        groups[chrom].append(L + "\n")
    i = 1
    for chrom, seqs in groups.items():
        out = Path(output_dir) / f"{prefix}#1#chr{i}.gff"
        out.write_text("".join(seqs))
        i += 1

def process_output_gff_files():
    for root, _, files in os.walk(source_dir):
        for fn in files:
            if fn.endswith(('.gff','.gff3')):
                fp = Path(root) / fn
                split_and_rename_gff(str(fp), splitgff_dir)

def split_gff_by_c_blocks(gff_path, output_dir):
    """按标准染色体（chr+数字）拆分GFF文件"""
    chr_groups = defaultdict(list)
    
    for line in Path(gff_path).read_text().splitlines():
        if line.startswith("#"):
            continue
        
        chrom = line.split('\t')[0]
        # 仅处理chr开头且带数字的染色体（如chr1, chr02）
        if re.match(r'chr\d+', chrom):
            chr_groups[chrom].append(line + "\n")
    
    # 按染色体名称排序后输出
    for idx, (chrom, lines) in enumerate(
        sorted(chr_groups.items(), key=lambda x: int(re.search(r'\d+', x[0]).group())), 
        start=1
    ):
        bed_dir = Path(output_dir) / f"{idx}.bed"
        bed_dir.mkdir(parents=True, exist_ok=True)
        (bed_dir / "yuangene.gff3").write_text("".join(lines))

def split_fna_by_chromosomes(fna_path, output_dir):
    """同样只处理标准染色体的FNA文件"""
    chrom_records = {}
    for record in SeqIO.parse(fna_path, "fasta"):
        if re.match(r'chr\d+', record.id):
            chrom_num = int(re.search(r'\d+', record.id).group())
            chrom_records[chrom_num] = record
    
    # 按染色体编号排序输出
    for chrom_num in sorted(chrom_records.keys()):
        bed_dir = Path(output_dir) / f"{chrom_num}.bed"
        bed_dir.mkdir(parents=True, exist_ok=True)
        SeqIO.write(chrom_records[chrom_num], bed_dir / "yuangene.fa", "fasta")
def generate_gffindex():
    print("🔍 生成 gffindex...")
    for fa in glob.glob(fa_pattern):
        fn   = Path(fa).stem
        outp = gffindex_dir / f"{fn}.txt"
        with open(fa) as fin, open(outp,'w') as fout:
            for L in fin:
                if L.startswith('>'):
                    fout.write(L)
    print("✅ gffindex 生成完成")

# ========== 新增注释处理功能 ==========
def load_annotations(index_file):
    """加载注释文件并返回一个字典，键为前缀，值为注释行的列表"""
    annotations = {}
    with open(index_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # 提取第一个#之前的内容作为前缀
                prefix = line.split('#')[0][1:]
                if prefix not in annotations:
                    annotations[prefix] = []
                annotations[prefix].append(line.strip())
    return annotations

def process_sp_file(sp_file, output_file, annotations):
    """处理单个 sp 文件，添加注释并输出到新文件"""
    prefix = sp_file.stem
    
    if prefix not in annotations:
        print(f"⚠️ 跳过 {sp_file}, 无对应注释")
        return
    
    annotation_lines = annotations[prefix]
    annotation_index = 0
    prev_start_pos_length = None
    
    with open(sp_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
                
            # 提取坐标信息
            chr_id = parts[0]
            node_id = parts[1]
            start_pos = parts[2]
            end_pos = parts[3]
            
            # 检查起始位置字符串长度
            current_start_pos_length = len(start_pos)
            
            # 当坐标长度变小时切换到下一个注释
            if prev_start_pos_length is not None and current_start_pos_length < prev_start_pos_length:
                annotation_index += 1
                if annotation_index >= len(annotation_lines):
                    print(f"⚠️ {sp_file} 的注释行不足，已使用最后一个注释")
                    annotation_index = len(annotation_lines) - 1
            
            # 更新前一个位置长度
            prev_start_pos_length = current_start_pos_length
            
            # 获取当前注释
            current_annotation = annotation_lines[annotation_index]
            
            # 写入新行（添加注释到第5列）
            outfile.write(f"{chr_id}\t{node_id}\t{start_pos}\t{end_pos}\t{current_annotation}\n")

def process_sp_files():
    """处理所有sp目录中的文件并添加注释"""
    print("\n🔍 开始处理注释添加...")
    
    for i in range(100):  # 假设最多100个子目录
        # 获取索引文件
        index_file = gffindex_dir / f"{i}all_genomes.txt"
        if not index_file.exists():
            continue
        
        # 加载注释
        annotations = load_annotations(index_file)
        
        # 查找对应的sp目录
        sp_dir = output_base / f"{i}.bed" / f"sp{i}"
        if not sp_dir.exists():
            continue
            
        # 处理每个bed文件
        for bed_file in sp_dir.glob("*.bed"):
            # 创建输出目录
            out_dir = annotated_dir / f"{i}.bed"
            out_dir.mkdir(exist_ok=True)
            
            # 处理文件
            process_sp_file(bed_file, out_dir / bed_file.name, annotations)
            
            print(f"✅ 已处理: {bed_file} → {out_dir / bed_file.name}")
    
    print("🎉 注释添加完成")

# ========== 主流程 ==========
def main():
    # 原有处理流程
    process_output_gff_files()
    split_gff_by_c_blocks(input_gff, output_base)
    split_fna_by_chromosomes(input_fna, output_base)
    unify_chr_numbers()
    generate_gffindex()
    
    # 新增注释处理
    process_sp_files()
    
    print("\n🎉 全部处理完成！")

if __name__ == "__main__":
    main()