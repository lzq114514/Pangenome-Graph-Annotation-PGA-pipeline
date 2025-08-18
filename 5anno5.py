#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
import time
import glob
import shutil
import subprocess
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from Bio import SeqIO

# ========= 用户可编辑的配置 =========
# SBATCH 参数配置 - 用户可根据需要修改这些值
SBATCH_JOB_NAME = "pep"                   # 作业名称
SBATCH_PARTITION = "hebhcnormal01"        # 使用的分区
SBATCH_NODES = 1                          # 节点数
SBATCH_NTASKS_PER_NODE = 60               # 每个节点的任务数
SBATCH_THREADS = 16                       # miniprot 使用的线程数

# ========= 与第一个脚本相同的输入接口 =========
if len(sys.argv) != 2:
    print("用法: python script.py <输入基础目录>")
    print("示例: python script.py /public/home/.../malingshufangfa2")
    sys.exit(1)

input_root   = sys.argv[1].rstrip(os.sep)
workflow2_dir = os.path.join(input_root+'workflow2')
workflow3_dir = os.path.join(input_root+'workflow3')

# workflow3 子目录
sp18_root    = os.path.join(workflow3_dir, "sp18")
getquery_dir = os.path.join(workflow3_dir, "getquery")
splitgff_dir = os.path.join(workflow3_dir, "splitgff")
qujian_dir   = os.path.join(workflow3_dir, "qujian")

for d in (getquery_dir, qujian_dir):
    os.makedirs(d, exist_ok=True)

MAX_PROCS = min(cpu_count(), 4)

PREFIX = "shuidao1#1#"
# ========== 原版GFF修复函数（放在最前面） ==========
def fix_gff(input_file, output_file):
    """
    完全按照您提供的原版修复函数
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # 写入 GFF3 版本声明
        print("##gff-version 3", file=f_out)
        
        gene_count = 1
        exon_count = 1           # 用于生成 exon ID
        current_parent = None    # 当前 mRNA ID，用于后续 CDS/exon/stop_codon
        
        for line in f_in:
            # 跳过注释、空行和 PAF 头部
            if line.startswith("#") or not line.strip() or line.startswith("##PAF"):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            seqid, source, feature, start, end, score, strand, phase, attr = fields
            
            # ----- mRNA -----
            if feature == "mRNA":
                exon_count = 1  # 为新的 transcript 重置 exon 计数
                
                # 提取或生成 mRNA ID
                m = re.search(r'ID=([^;\s]+)', attr)
                if m:
                    mrna_id = m.group(1)
                else:
                    mrna_id = f"MP{gene_count:06d}"
                    attr = f"ID={mrna_id};" + attr
                
                # 输出 gene 行
                gene_id = f"gene-{mrna_id}"
                gene_attr = f"ID={gene_id};Name={gene_id}"
                print("\t".join([seqid, source, "gene", start, end, score, strand, ".", gene_attr]), file=f_out)
                
                # 更新并输出 mRNA 行，确保有 Parent
                attr = re.sub(r'ID=[^;\s]+', f"ID={mrna_id};Parent={gene_id}", attr)
                print("\t".join([seqid, source, "mRNA", start, end, score, strand, ".", attr]), file=f_out)
                
                current_parent = mrna_id
                gene_count += 1
                continue  # 直接进入下一行，不走下面的其他分支
            
            # ----- exon （保留输入中的 exon） -----
            if feature == "exon" and current_parent:
                # 如果没有 ID，就生成；并确保 Parent
                if 'ID=' not in attr:
                    exon_id = f"exon-{current_parent}-{exon_count}"
                    exon_count += 1
                    attr = f"ID={exon_id};Parent={current_parent}"
                else:
                    # 补全 Parent
                    if 'Parent=' in attr:
                        attr = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                    else:
                        attr += f";Parent={current_parent}"
                print("\t".join([seqid, source, "exon", start, end, score, strand, phase, attr]), file=f_out)
                continue
            
            # ----- CDS （并自动生成 exon） -----
            if feature == "CDS" and current_parent:
                # 1) 自动生成 exon
                exon_id = f"exon-{current_parent}-{exon_count}"
                exon_count += 1
                exon_attr = f"ID={exon_id};Parent={current_parent}"
                print("\t".join([seqid, source, "exon", start, end, ".", strand, ".", exon_attr]), file=f_out)
                
                # 2) 输出 CDS
                if 'Parent=' in attr:
                    attr_cds = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                else:
                    attr_cds = f"Parent={current_parent};" + attr
                print("\t".join([seqid, source, "CDS", start, end, score, strand, phase, attr_cds]), file=f_out)
                continue
            
            # ----- stop_codon -----
            if feature == "stop_codon" and current_parent:
                if 'Parent=' in attr:
                    attr_stop = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                else:
                    attr_stop = f"Parent={current_parent};" + attr
                print("\t".join([seqid, source, "stop_codon", start, end, score, strand, phase, attr_stop]), file=f_out)
                continue

# ========= STEP1: BED -> getquery =========
def process_bed_file(mod2_dir, sp18_dir, out_dir):
    print(f"[STEP1] 处理 {mod2_dir}")
    os.makedirs(out_dir, exist_ok=True)
    for bed in glob.glob(os.path.join(mod2_dir, "*.bed")):
        start = time.time()
        base = os.path.splitext(os.path.basename(bed))[0]
        sp18_file = os.path.join(sp18_dir, re.sub(r"_hap1$", "", base) + ".bed")
        if not os.path.exists(sp18_file):
            print(f"  跳过，未找到 {sp18_file}")
            continue

        # 读取 sp18
        locs = defaultdict(dict)
        with open(sp18_file) as fh:
            for L in fh:
                p = L.split()
                if len(p) < 5: continue
                cid, nid, st, ed, *rest = p
                locs[cid][nid] = (st, ed, rest)

        out_file = os.path.join(out_dir, os.path.basename(bed))
        with open(bed) as fin, open(out_file, 'w') as fout:
            for L in fin:
                p = L.split()
                if len(p) < 2: continue
                cid, *nids = p
                if cid not in locs: continue
                for nid in nids:
                    if nid in locs[cid]:
                        st, ed, rest = locs[cid][nid]
                        fout.write(f"{cid}\t{st}\t{ed}\t{' '.join(rest)}\n")

        print(f"  完成 {os.path.basename(bed)}，用时 {time.time()-start:.2f}s")

# ========= STEP2: getquery -> qujian =========
def process_qujian_folder(num):
    bed_folder = os.path.join(getquery_dir, f"{num}.bed")
    out_folder = os.path.join(qujian_dir, num)
    os.makedirs(out_folder, exist_ok=True)
    print("\n" + "="*40)
    print(f"[STEP2] 处理 getquery -> qujian: {num}")
    print(f"  输入 BED 目录: {bed_folder}")
    print(f"  输出 GFF 目录: {out_folder}")

    # 1. 加载并合并区间
    bed_mapping = defaultdict(list)
    for bed_file in glob.glob(os.path.join(bed_folder, "*.bed")):
        for line in open(bed_file):
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            chrom = parts[3].split(">")[-1].replace("chr0", "chr")
            bed_mapping[chrom].append((start, end))

    print(f"  载入染色体: {sorted(bed_mapping.keys())}")

    # 合并区间
    optimized = {}
    for chrom, intervals in bed_mapping.items():
        intervals.sort(key=lambda x: x[0])
        merged = []
        for cur in intervals:
            if not merged or cur[0] > merged[-1][1]:
                merged.append([cur[0], cur[1]])
            else:
                merged[-1][1] = max(merged[-1][1], cur[1])
        optimized[chrom] = merged
        print(f"  {chrom.ljust(6)}: 合并为 {len(merged)} 个区间")

    # 2. 扫描 GFF 并输出
    processed = 0
    if not os.path.isdir(splitgff_dir):
        print(f"  警告: 找不到 splitgff 目录: {splitgff_dir}，跳过此 num 的 GFF 扫描。")
    else:
        for gff_file in sorted(os.listdir(splitgff_dir)):
            if not (gff_file.endswith(".gff") or gff_file.endswith(".gff3")):
                continue
            chrom = os.path.splitext(gff_file)[0].replace("chr0", "chr")
            if chrom not in optimized:
                continue

            ivs = optimized[chrom]
            in_path  = os.path.join(splitgff_dir, gff_file)
            out_path = os.path.join(out_folder, gff_file)
            with open(in_path) as fin, open(out_path, 'w') as fout:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    cols = line.split("\t")
                    try:
                        gff_s, gff_e = int(cols[3]), int(cols[4])
                    except:
                        continue

                    # 二分查找判断是否落入任一区间
                    lo, hi, hit = 0, len(ivs)-1, False
                    while lo <= hi:
                        mid = (lo + hi) // 2
                        b0, b1 = ivs[mid]
                        if gff_e < b0:
                            hi = mid - 1
                        elif gff_s > b1:
                            lo = mid + 1
                        else:
                            hit = True
                            break
                    if hit:
                        fout.write(line)

            processed += 1
            print(f"  输出: {gff_file} → {out_path}")

    print(f"[STEP2] 完成 {num}，共处理 {processed} 个 GFF 文件")
    print("="*40)

# ========= workflow3 -> workflow4 helper functions =========
def detect_workflow4_dir(workflow3_dir):
    if not workflow3_dir:
        raise ValueError("请提供 workflow3 路径。")
    workflow3_dir = os.path.abspath(workflow3_dir)
    low = workflow3_dir.lower()
    idx = low.rfind('workflow3')
    if idx != -1:
        prefix = workflow3_dir[:idx]
        wf4 = prefix + 'workflow4'
    else:
        wf4 = os.path.join(os.path.dirname(workflow3_dir), 'workflow4')
    os.makedirs(wf4, exist_ok=True)
    return wf4

def detect_workflow1_dir(workflow3_dir):
    workflow3_dir = os.path.abspath(workflow3_dir)
    low = workflow3_dir.lower()
    idx = low.rfind('workflow3')
    if idx != -1:
        prefix = workflow3_dir[:idx]
        wf1 = prefix + 'workflow1'
        if os.path.isdir(wf1):
            return wf1
    parent_wf1 = os.path.join(os.path.dirname(workflow3_dir), 'workflow1')
    if os.path.isdir(parent_wf1):
        return parent_wf1
    return None

# ========= 合并 qujian 下按前缀的 GFF =========
def merge_gff_by_prefix(qujian_dir, output_dir):
    if not os.path.isdir(qujian_dir):
        print(f"警告: 找不到 qujian 目录: {qujian_dir}，跳过合并。")
        return []

    subdirs = [d for d in os.listdir(qujian_dir) if os.path.isdir(os.path.join(qujian_dir, d))]
    prefix_map = defaultdict(list)
    for subdir in subdirs:
        subdir_path = os.path.join(qujian_dir, subdir)
        gff_files = glob.glob(os.path.join(subdir_path, "*.gff*"))
        for gff_file in gff_files:
            filename = os.path.basename(gff_file)
            if '#' in filename:
                prefix = filename.split('#', 1)[0]
            else:
                prefix = os.path.splitext(filename)[0]
            prefix_map[prefix].append(gff_file)

    if not prefix_map:
        print(f"未在 {qujian_dir} 下发现任何 GFF 文件。")
        return []

    written_files = []
    for prefix, file_list in prefix_map.items():
        file_list.sort()
        output_file = os.path.join(output_dir, f"{prefix}.gff")
        with open(output_file, "w") as out_f:
            if file_list:
                with open(file_list[0]) as first_f:
                    out_f.write(first_f.read())
            for file_path in file_list[1:]:
                with open(file_path) as in_f:
                    for line in in_f:
                        if not line.startswith("#"):
                            out_f.write(line)
        written_files.append(output_file)
        print(f"合并完成: {len(file_list)} 个文件 -> {output_file}")
    return written_files

# ========= 复制 .fna =========
def copy_fna_from_workflow1(workflow1_dir, workflow4_dir, overwrite=False):
    if not workflow1_dir or not os.path.isdir(workflow1_dir):
        print("未找到 workflow1 目录，跳过 .fna 复制。")
        return []

    os.makedirs(workflow4_dir, exist_ok=True)

    fna_patterns = [os.path.join(workflow1_dir, '**', '*.fna'),
                    os.path.join(workflow1_dir, '**', '*.fa'),
                    os.path.join(workflow1_dir, '**', '*.fasta')]
    copied = []
    seen_names = set()
    for pat in fna_patterns:
        for src in glob.glob(pat, recursive=True):
            base = os.path.basename(src)
            if base in seen_names:
                print(f"跳过重复文件名: {src}")
                continue
            seen_names.add(base)
            dst = os.path.join(workflow4_dir, base)
            if os.path.exists(dst) and not overwrite:
                print(f"目标已存在，跳过: {dst}")
                continue
            shutil.copy2(src, dst)
            copied.append(dst)
            print(f"已复制 .fna: {src} -> {dst}")
    return copied

# ========= 补齐 workflow1 中缺失但 workflow4 中不存在的 gff =========
def copy_missing_gffs_from_workflow1(workflow1_dir, workflow4_dir, overwrite=False):
    if not workflow1_dir or not os.path.isdir(workflow1_dir):
        print("未找到 workflow1 目录，跳过 gff 补齐。")
        return []

    existing = set()
    for f in glob.glob(os.path.join(workflow4_dir, '*.gff*')):
        existing.add(os.path.basename(f))

    copied = []
    for src in glob.glob(os.path.join(workflow1_dir, '**', '*.gff*'), recursive=True):
        name = os.path.basename(src)
        dst = os.path.join(workflow4_dir, name)
        if name in existing and os.path.exists(dst) and not overwrite:
            continue
        shutil.copy2(src, dst)
        copied.append(dst)
        print(f"已复制 gff: {src} -> {dst}")
    return copied

# ========= 提取蛋白并合并 =========
def extract_and_merge_proteins(input_dir):
    """
    1. 在 input_dir 下查找成对的 fna + gff（基于前缀）
    2. 使用 gffread 提取蛋白序列
    3. 合并所有蛋白到 input_dir/miniprotzhushi.fa
    """
    output_dir = os.path.join(input_dir, "protein_results")
    os.makedirs(output_dir, exist_ok=True)

    file_groups = defaultdict(dict)

    # 寻找 fna (*.fna, *.fa, *.fasta) 和 gff 文件
    for fna_ext in ("*.fna", "*.fa", "*.fasta"):
        for fna_file in glob.glob(os.path.join(input_dir, fna_ext)):
            base = os.path.basename(fna_file)
            prefix = base.split(".")[0]
            file_groups[prefix]["fna"] = fna_file

    for gff_file in glob.glob(os.path.join(input_dir, "*.gff*")):
        base = os.path.basename(gff_file)
        prefix = base.split(".")[0]
        file_groups[prefix]["gff"] = gff_file

    all_proteins = []
    for prefix, files in file_groups.items():
        if "fna" not in files or "gff" not in files:
            print(f"跳过 {prefix} - 缺少 fna 或 gff 文件")
            continue

        fna_path = files["fna"]
        gff_path = files["gff"]
        protein_file = os.path.join(output_dir, f"{prefix}_protein.fa")

        print(f"处理: {prefix}")
        print(f"  FNA: {fna_path}")
        print(f"  GFF: {gff_path}")

        try:
            subprocess.run([
                "gffread",
                "-y", protein_file,
                "-g", fna_path,
                gff_path
            ], check=True)
            all_proteins.append(protein_file)
            print(f"  成功生成: {protein_file}")
        except subprocess.CalledProcessError as e:
            print(f"  gffread 失败: {e}")
            continue
        except FileNotFoundError:
            print("  错误: gffread 未找到。请确保 gffread 已安装且在 PATH 中。")
            return

    # 合并
    if all_proteins:
        cat_file = os.path.join(input_dir, "miniprotzhushi.fa")
        print(f"\n合并所有蛋白质序列到: {cat_file}")
        with open(cat_file, "w") as out_handle:
            for protein_file in sorted(all_proteins):
                for record in SeqIO.parse(protein_file, "fasta"):
                    SeqIO.write(record, out_handle, "fasta")
        print(f"合并完成! 共 {len(all_proteins)} 个蛋白质文件")
    else:
        print("没有生成任何蛋白质文件")

def run_miniprot_and_fix(workflow4_dir):
    """运行miniprot并修复GFF"""
    # 1. 准备文件路径
    genome_fasta = os.path.join(workflow4_dir, "shuidao1.fna")
    index_file = os.path.join(workflow4_dir, "genome.mpi")
    protein_fasta = os.path.join(workflow4_dir, "miniprotzhushi.fa")
    raw_gff = os.path.join(workflow4_dir, "aln.gff")
    fixed_gff = os.path.join(workflow4_dir, "pga_anno.gff")

    # 2. 检查必要文件
    if not os.path.exists(genome_fasta):
        print(f"错误: 未找到基因组文件 {genome_fasta}")
        return False
    if not os.path.exists(protein_fasta):
        print(f"错误: 未找到蛋白文件 {protein_fasta}")
        return False

    # 3. 生成sbatch脚本 - 使用用户配置的变量
    sbatch_content = f"""#!/bin/bash
#SBATCH --job-name={SBATCH_JOB_NAME}
#SBATCH --partition={SBATCH_PARTITION}
#SBATCH --nodes={SBATCH_NODES}
#SBATCH --ntasks-per-node={SBATCH_NTASKS_PER_NODE}
#SBATCH --error=%j.err
#SBATCH --output=%j.out

# 生成索引
miniprot -t16 -d {index_file} {genome_fasta}

# 运行比对
miniprot -t16 --gff {index_file} {protein_fasta} > {raw_gff}
"""
    sbatch_path = os.path.join(workflow4_dir, "run_miniprot.sbatch")
    with open(sbatch_path, 'w') as f:
        f.write(sbatch_content)

    # 4. 提交作业
    try:
        print("提交miniprot作业...")
        result = subprocess.run(
            ['sbatch', sbatch_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        job_id = result.stdout.split()[-1]
        print(f"作业已提交，ID: {job_id}")

        # 5. 等待作业完成（简化版）
        print("等待miniprot完成（约5分钟）...")
        time.sleep(300)

        # 6. 修复GFF
        if os.path.exists(raw_gff):
            print("开始修复GFF文件...")
            fix_gff(raw_gff, fixed_gff)  # 现在可以正确找到fix_gff函数
            print(f"GFF修复完成: {fixed_gff}")
            return True
        else:
            print(f"错误: 未生成原始GFF文件 {raw_gff}")
            return False

    except subprocess.CalledProcessError as e:
        print(f"作业提交失败: {e.stderr}")
        return False

def main():
    # 初始化路径
    if len(sys.argv) != 2:
        print("用法: python script.py <输入基础目录>")
        sys.exit(1)

    input_root = sys.argv[1].rstrip(os.sep)
    workflow2_dir = os.path.join(input_root+'workflow2')
    workflow3_dir = os.path.join(input_root+'workflow3')
    sp18_root = os.path.join(workflow3_dir, "sp18")
    getquery_dir = os.path.join(workflow3_dir, "getquery")
    splitgff_dir = os.path.join(workflow3_dir, "splitgff")
    qujian_dir = os.path.join(workflow3_dir, "qujian")

    # 创建必要目录
    for d in (getquery_dir, qujian_dir):
        os.makedirs(d, exist_ok=True)

    # 1. 处理BED文件
    tasks = []
    for bd in glob.glob(os.path.join(workflow3_dir, "*.bed")):
        if os.path.isdir(bd):
            m2 = os.path.join(bd, "modify2")
            s18 = os.path.join(sp18_root, os.path.basename(bd))
            out = os.path.join(getquery_dir, os.path.basename(bd))
            if os.path.isdir(m2) and os.path.isdir(s18):
                tasks.append((m2, s18, out))
    
    if tasks:
        with Pool(MAX_PROCS) as p:
            p.starmap(process_bed_file, tasks)

    # 2. 处理区间文件夹
    nums = sorted([os.path.splitext(os.path.basename(d))[0]
                  for d in glob.glob(os.path.join(getquery_dir, "*.bed"))
                  if os.path.isdir(d)], key=int)
    if nums:
        with Pool(MAX_PROCS) as p:
            p.map(process_qujian_folder, nums)

    # 3. 创建workflow4目录
    workflow4_dir = detect_workflow4_dir(workflow3_dir)
    os.makedirs(workflow4_dir, exist_ok=True)
    print(f"输出目录: {workflow4_dir}")

    # 4. 合并GFF文件
    merge_gff_by_prefix(qujian_dir, workflow4_dir)

    # 5. 复制FNA文件
    workflow1_dir = os.path.join(input_root+'workflow1')
    if os.path.exists(workflow1_dir):
        copy_fna_from_workflow1(workflow1_dir, workflow4_dir)

    # 6. 补齐workflow1中的gff文件到workflow4 (添加这一步)
    if os.path.exists(workflow1_dir):
        print("\n[补齐步骤] 复制workflow1中的gff文件到workflow4")
        copy_missing_gffs_from_workflow1(workflow1_dir, workflow4_dir)
    else:
        print("\n[补齐步骤] 未找到workflow1目录，跳过gff文件补齐")

    # 7. 提取合并蛋白质
    extract_and_merge_proteins(workflow4_dir)

    # 8. 运行miniprot并修复GFF
    if not run_miniprot_and_fix(workflow4_dir):
        print("miniprot流程执行失败")
    else:
        print("所有流程成功完成")

if __name__ == "__main__":
    main()