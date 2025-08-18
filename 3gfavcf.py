#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import glob
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd

# ========== 参数解析 ==========
if len(sys.argv) < 2:
    print("用法: python script.py <输入目录> [GFF_FNA_前缀] [GFF标签前缀]")
    sys.exit(1)

input_root     = sys.argv[1]
gff_fna_prefix = sys.argv[2] if len(sys.argv) > 2 else None
gff_tag_prefix = sys.argv[3] if len(sys.argv) > 3 else "NC"

# ========== 确保 workflow3 输出目录存在 ==========
workflow3_dir = input_root + "workflow3"
os.makedirs(workflow3_dir, exist_ok=True)


# ========== 第一步：GFF ➝ BED ==========
print("Step I: GFF ➝ BED 转换…")
input_gff_dir  = os.path.join(input_root, "output")
gff_output_dir = os.path.join(workflow3_dir, "gff1")
os.makedirs(gff_output_dir, exist_ok=True)

for gff_file in glob.glob(os.path.join(input_gff_dir, "*.gff")):
    out_bed = os.path.join(gff_output_dir, os.path.basename(gff_file).replace(".gff", ".bed"))
    with open(gff_file, 'r') as fin, open(out_bed, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                chrom = parts[0]
                start = int(parts[3]) - 1
                end   = parts[4]
                gene  = parts[8].split(';')[0].replace("ID=", "")
                fout.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
    print(f"  → {os.path.basename(gff_file)} 转换为 {os.path.basename(out_bed)}")

# ========== 第二步：VCF + BED 处理 ==========
print("Step II: VCF + BED 处理…")
vcf_pattern = input_root + "workflow2" + "/*.pggb_out/chr*.vcf"
vcf_files   = glob.glob(vcf_pattern)
if not vcf_files:
    print("⚠️ 未找到任何 VCF 文件，退出")
    sys.exit(1)

# 取第一个 GFF ➝ BED 作为参考
ref_bed_files = glob.glob(os.path.join(gff_output_dir, "*.bed"))
if not ref_bed_files:
    print("⚠️ 未找到参考 BED，退出")
    sys.exit(1)
ref_bed = ref_bed_files[0]
temp_ref = os.path.join(workflow3_dir, "temp_ref.bed")
shutil.copy(ref_bed, temp_ref)

processed_dirs = 0
for vcf in vcf_files:
    m = re.search(r'chr(\d+)\.vcf$', vcf)
    if not m:
        continue
    idx   = m.group(1)
    outd  = os.path.join(workflow3_dir, f"{idx}.bed")
    os.makedirs(outd, exist_ok=True)
    qbed  = os.path.join(outd, "query.bed")
    try:
        subprocess.run(["perl", "vcfbed.pl", vcf, temp_ref, qbed], check=True)
        print(f"  ✔ 处理 {os.path.basename(vcf)} -> query.bed")
        processed_dirs += 1
    except subprocess.CalledProcessError as e:
        print(f"  ❌ 处理失败 {vcf}: {e}")
# 删除临时参考
os.remove(temp_ref)

# ========== 第三步：query.bed ➝ *_hap1.bed ==========
print("Step III: 拆分 hap1…")
for i in range(1, 100):
    d = os.path.join(workflow3_dir, f"{i}.bed")
    qp = os.path.join(d, "query.bed")
    if not os.path.isfile(qp):
        continue
    groups = defaultdict(list)
    with open(qp) as f:
        for L in f:
            if L.strip():
                chrom = L.split('\t',1)[0]
                groups[f"{chrom}_hap1"].append(L)
    for name, lines in groups.items():
        with open(os.path.join(d, f"{name}.bed"), 'w') as fo:
            fo.writelines(lines)

# ========== 第四步：GFA ➝ BED ==========
print("Step IV: GFA ➝ BED…")
gfa_files = glob.glob(input_root + "workflow2" + "/*.pggb_out/*.gfa")
for gfa in gfa_files:
    pid = re.search(r"(\d+)\.pggb_out", os.path.dirname(gfa))
    if not pid:
        continue
    outd = os.path.join(workflow3_dir, f"{pid.group(1)}.bed")
    os.makedirs(outd, exist_ok=True)
    outf = os.path.join(outd, os.path.basename(gfa).replace(".gfa", ".bed"))
    with open(outf, 'w') as fo:
        subprocess.run(
            ["perl", "/public/home/acfurbn1nz/pan/gfa2bed.pl", gfa],
            check=True, stdout=fo
        )

# ========== 第五步：重命名 .final.bed ==========
print("Step V: 重命名 .final.bed…")
for d in glob.glob(os.path.join(workflow3_dir, "*.bed")):
    if os.path.isdir(d):
        cnt = 0
        for f in sorted(glob.glob(os.path.join(d, "*.final.bed"))):
            os.rename(f, os.path.join(d, f"chr{cnt}.bed"))
            cnt += 1

# ========== 第六步：清理 hap*.bed ➝ modify2 ==========
print("Step VI: 清理 hap*.bed…")
for d in glob.glob(os.path.join(workflow3_dir, "*.bed")):
    if os.path.isdir(d):
        md = os.path.join(d, "modify2")
        os.makedirs(md, exist_ok=True)
        for f in os.listdir(d):
            if 'hap' in f and f.endswith('.bed'):
                with open(os.path.join(d, f)) as fin, open(os.path.join(md, f), 'w') as fout:
                    for L in fin:
                        fout.write(re.sub(r'[<>]', ' ', L.strip()) + "\n")

# ========== 第七步：处理 chr0.bed 的第一列去掉 '#' 及其后内容 ==========
print("Step VII: 清理 chr0.bed 首列 ‘#’ 标签…")
for i in range(100):
    bd = os.path.join(workflow3_dir, f"{i}.bed")
    f0 = os.path.join(bd, "chr0.bed")
    if os.path.isfile(f0):
        lines = open(f0).read().splitlines()
        new = []
        for L in lines:
            cols = re.split(r'\s+', L)
            cols[0] = re.sub(r'#.*','', cols[0])
            new.append("\t".join(cols))
        open(f0, 'w').write("\n".join(new) + "\n")

# ========== 第八步：GFF 分组 & FNA 分组 ==========
print("Step VIII: GFF & FNA 分组…")
def split_gff_by_prefix(gff_path, out_dir, prefix):
    print(f"  🔍 处理GFF文件: {gff_path}")
    print(f"  使用前缀: '{prefix}'")
    groups, curr, cid = [], [], None
    total_lines = 0
    processed_lines = 0
    
    try:
        with open(gff_path, 'r') as f:
            for line_num, L in enumerate(f, 1):
                total_lines += 1
                if L.startswith("#"):
                    continue
                chrom = L.split('\t', 1)[0]
                if chrom.upper().startswith(prefix.upper()):
                    processed_lines += 1
                    if chrom != cid:
                        if curr: 
                            groups.append(curr)
                        curr, cid = [L], chrom
                    else:
                        curr.append(L)
                elif curr:
                    curr.append(L)
        if curr: 
            groups.append(curr)
    except Exception as e:
        print(f"  ❌ 处理GFF文件时出错: {e}")
        return
    
    print(f"  总行数: {total_lines}, 处理行数: {processed_lines}")
    print(f"  找到 {len(groups)} 个分组")
    
    if not groups:
        print("  ⚠️ 警告: 未找到任何匹配的分组")
        return
    
    for idx, grp in enumerate(groups, 1):
        gd = Path(out_dir) / f"{idx}.bed"
        gd.mkdir(parents=True, exist_ok=True)
        output_file = gd / "yuangene.gff3"
        with open(output_file, 'w') as fo:
            fo.write("##gff-version 3\n")
            fo.writelines(grp)
        print(f"  分组 {idx}: 包含 {len(grp)} 行 → {output_file}")

def split_fna_by_chr(fna_path, out_dir):
    print(f"  🔍 处理FNA文件: {fna_path}")
    try:
        content = open(fna_path).read().strip()
        entries = [e for e in content.split('>') if e]
        print(f"  找到 {len(entries)} 个FASTA条目")
        
        cmap = {}
        for ent in entries:
            lines = ent.split('\n')
            hdr = lines[0]
            seq = "\n".join(lines[1:])
            
            # 调试: 打印原始标题
            print(f"    处理条目: {hdr[:50]}{'...' if len(hdr) > 50 else ''}")
            
            # 修改：包括对 X 和 Y 的处理
            m = re.search(r'chromosome\s*(\d+|X|Y)', hdr, re.IGNORECASE)
            if not m:
                m = re.search(r'chr\s*(\d+|X|Y)', hdr, re.IGNORECASE)
            if not m:
                m = re.search(r'(\d+|X|Y)$', hdr)
                
            if m:
                n = m.group(1)
                cmap.setdefault(n, []).append(f">{hdr}\n{seq}")
                print(f"      匹配染色体: {n}")
            else:
                print(f"    ⚠️ 无法识别染色体编号: {hdr[:50]}{'...' if len(hdr) > 50 else ''}")
        
        if not cmap:
            print("  ❌ 错误: 未找到任何有效的染色体条目")
            return
        
        print(f"  找到 {len(cmap)} 个染色体分组")
        for n, recs in cmap.items():
            od = Path(out_dir) / f"{n}.bed"
            od.mkdir(parents=True, exist_ok=True)
            output_file = od / "yuangene.fa"
            with open(output_file, 'w') as f:
                f.write("\n".join(recs))
            print(f"    染色体 {n}: {len(recs)} 条序列 → {output_file}")
    except Exception as e:
        print(f"  ❌ 处理FNA文件时出错: {e}")


if gff_fna_prefix:
    # 修正路径构建 - 文件在输入目录的 malingshu2 子目录中
    malingshu2_dir =  input_root + "workflow2"
    
    # 检查 malingshu2 目录是否存在
    if not os.path.isdir(malingshu2_dir):
        print(f"  ❌ 错误: malingshu2 目录不存在: {malingshu2_dir}")
        print(f"  输入目录内容: {os.listdir(input_root)}")
    else:
        gff_path = os.path.join(input_root + "workflow1", f"{gff_fna_prefix}.gff")
        fna_path = os.path.join(input_root + "workflow1", f"{gff_fna_prefix}.fna")
        
        print(f"  GFF路径: {gff_path}")
        print(f"  FNA路径: {fna_path}")
        
        if os.path.isfile(gff_path):
            print("  ✅ GFF文件存在")
            split_gff_by_prefix(gff_path, workflow3_dir, gff_tag_prefix)
        else:
            print(f"  ❌ GFF文件不存在: {gff_path}")
            # 列出输入目录内容以帮助调试
            print(f"  malingshu2 目录内容: {os.listdir(malingshu2_dir)}")
        
        if os.path.isfile(fna_path):
            print("  ✅ FNA文件存在")
            split_fna_by_chr(fna_path, workflow3_dir)
        else:
            print(f"  ❌ FNA文件不存在: {fna_path}")
            # 列出输入目录内容以帮助调试
            print(f"  malingshu2 目录内容: {os.listdir(malingshu2_dir)}")
else:
    print("  ⚠️ 未提供GFF/FNA前缀，跳过第八步")
# ========== 第九步：按 chr0.bed 首列 ID 拆分 ==========
print("Step IX: 拆分 chr0.bed by ID…")
for bd in os.listdir(workflow3_dir):
    full = os.path.join(workflow3_dir, bd)
    if os.path.isdir(full) and bd.endswith(".bed"):
        idx = bd.replace(".bed","")
        outd = os.path.join(full, f"sp{idx}")
        os.makedirs(outd, exist_ok=True)
        f0 = os.path.join(full, "chr0.bed")
        if not os.path.isfile(f0):
            print(f"  ⚠️ 跳过 {bd}：未找到 chr0.bed")
            continue
        df = pd.read_csv(f0, sep="\t", header=None, dtype=str)
        for gid, grp in df.groupby(0):
            grp.to_csv(os.path.join(outd, f"{gid}.bed"), sep="\t", index=False, header=False)
        print(f"  ✔ 拆分完成: {bd} → sp{idx}/")

print(f"🎉 全部完成，共处理 0 个目录。")