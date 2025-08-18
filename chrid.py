#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import shutil
from pathlib import Path
from collections import defaultdict
import pandas as pd

# ===================== 用户配置区域 =====================
# EMAPPER_PATH = "/path/to/emapper.py"
# ======================================================

def extract_prefix(filename):
    """
    提取文件名前缀（连续的字母数字序列）。
    """
    match = re.match(r'^[A-Za-z0-9]+', filename)
    return match.group(0) if match else None


def process_annotation_pair(gff_path, anno_path, output_dir):
    """
    根据 annotation 文件过滤 GFF，输出到 output_dir。
    """
    with open(gff_path, 'r') as f:
        gff_lines = f.readlines()
    try:
        df = pd.read_csv(anno_path, sep='\t', header=None, skiprows=5)
    except pd.errors.ParserError:
        try:
            df = pd.read_csv(anno_path, sep=',', header=None, skiprows=5)
        except pd.errors.ParserError:
            print(f"无法解析文件: {anno_path}")
            return
    queries = df.iloc[:, 0].dropna().tolist()
    matched_ids = set()
    for line in gff_lines:
        id_match = re.search(r'ID=([^;]*)', line)
        if id_match:
            matched_ids.add(id_match.group(1))
        parent_match = re.search(r'Parent=([^;]*)', line)
        if parent_match and parent_match.group(1) in queries:
            matched_ids.add(parent_match.group(1))
    out_file = os.path.join(output_dir, os.path.basename(anno_path).replace('.annotations', '.gff'))
    with open(out_file, 'w') as out_f:
        for line in gff_lines:
            id_match = re.search(r'ID=([^;]*)', line)
            if id_match and id_match.group(1) in matched_ids:
                out_f.write(line)
                continue
            parent_match = re.search(r'Parent=([^;]*)', line)
            if parent_match and parent_match.group(1) in matched_ids:
                out_f.write(line)
    print(f"✅ 注释过滤完成: {out_file}")


def collect_files(input_dir, output_dir):
    """
    复制 .gff, .fna, .pep 文件到 workflow1 目录
    """
    for root, _, files in os.walk(input_dir):
        for fn in files:
            if fn.endswith(('.gff', '.fna')) or fn.endswith('.pep'):
                shutil.copy2(os.path.join(root, fn), os.path.join(output_dir, fn))
                print(f"已复制: {fn}")


def split_fna_by_chromosome(input_dir, output_dir):
    """
    根据 FASTA 标题中 'chromosome <label>' 后面的数字或 X/Y 分组。
    """
    pattern = re.compile(r'chromosome\s*(\d+|X|Y)\b', re.IGNORECASE)

    for f in Path(input_dir).glob('*.fna'):
        prefix = f.stem
        groups = defaultdict(list)
        lines = f.read_text().splitlines()
        header, seq_lines = None, []

        # 遍历每条记录
        for line in lines + ['>END']:
            if line.startswith('>'):
                if header:
                    m = pattern.search(header)
                    if m:
                        label = m.group(1).upper()
                        seq = ''.join(seq_lines)
                        groups[label].append((header, seq))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # 写文件
        for label, entries in groups.items():
            d = Path(output_dir) / label
            d.mkdir(parents=True, exist_ok=True)
            out_path = d / f"{prefix}.fna"
            with open(out_path, 'w') as w:
                for hdr, seq in entries:
                    w.write(f"{hdr}\n{seq}\n")
            print(f"✅ 写入 {label}/{prefix}.fna ({len(entries)} 条)")


def modify_fna_ids(split_dir):
    """
    对每个 prefix 维护独立的全局计数器，确保同一 prefix
    在所有子目录中的 .fna 文件里连续编号。
    """
    base = Path(split_dir)
    
    # === 修复：添加染色体排序逻辑 ===
    chr_dirs = []
    for d in base.iterdir():
        if d.is_dir():
            name = d.name.upper()
            if name.isdigit():
                chr_dirs.append((int(name), d))
            elif name in ["X", "Y"]:
                chr_dirs.append((100 if name == "X" else 101, d))
    chr_dirs.sort(key=lambda x: x[0])
    # =============================

    counters = defaultdict(int)

    for _, d in chr_dirs:
        for f in sorted(d.glob('*.fna')):
            prefix = extract_prefix(f.name)
            if not prefix:
                continue
            lines = f.read_text().splitlines(keepends=True)
            new_lines = []
            for ln in lines:
                if ln.startswith('>'):
                    counters[prefix] += 1
                    new_chr = str(counters[prefix]).zfill(2)
                    new_lines.append(f">{prefix}#1#chr{new_chr}\n")
                else:
                    new_lines.append(ln)
            f.write_text(''.join(new_lines))

    for prefix, count in counters.items():
        print(f"✅ {prefix} 最终计数到 chr{str(count).zfill(2)}")


def merge_chromosome_files(split_dir):
    """
    合并各染色体子目录下的 .fna 到 single-chromosome FASTA
    """
    base = Path(split_dir).parent
    
    # === 修复：添加染色体排序逻辑 ===
    chr_dirs = []
    for d in Path(split_dir).iterdir():
        if d.is_dir():
            name = d.name.upper()
            if name.isdigit():
                chr_dirs.append((int(name), d))
            elif name in ["X", "Y"]:
                chr_dirs.append((100 if name == "X" else 101, d))
    chr_dirs.sort(key=lambda x: x[0])
    # =============================
    
    for i, (_, d) in enumerate(chr_dirs, 1):
        out = base / f"{i}all_genomes.fa"
        with open(out, 'w') as w:
            for f in sorted(d.glob('*.fna')):
                w.write(f.read_text())
        print(f"✅ 合并完成：{out}")


def main(input_dir):
    output_dir       = os.path.join(input_dir + 'workflow1')
    final_output_dir = os.path.join(input_dir, 'output')
    workflow2_dir    = os.path.join(input_dir + 'workflow2')
    split_dir        = os.path.join(workflow2_dir, 'chromosomes')

    for d in (output_dir, final_output_dir, workflow2_dir, split_dir):
        os.makedirs(d, exist_ok=True)

    print("\n=== 收集文件 ===")
    collect_files(input_dir, output_dir)

    print("\n=== 注释过滤 (.annotations -> .gff) ===")
    for fn in os.listdir(output_dir):
        if fn.endswith('.annotations'):
            prefix = extract_prefix(fn)
            gff_fn = prefix + '.gff'
            gff_path = os.path.join(output_dir, gff_fn)
            anno_path = os.path.join(output_dir, fn)
            if os.path.isfile(gff_path):
                process_annotation_pair(gff_path, anno_path, final_output_dir)

    print("\n=== 分割 fna ===")
    split_fna_by_chromosome(output_dir, split_dir)

    print("\n=== 修改 ID（跨染色体连续计数） ===")
    modify_fna_ids(split_dir)

    print("\n=== 合并染色体文件 ===")
    merge_chromosome_files(split_dir)

    print("\n🎉 全部完成！")
    print(f"- workflow1 结果: {output_dir}")
    print(f"- 注释过滤结果: {final_output_dir}")
    print(f"- split fna 目录: {split_dir}")
    print(f"- 合并后文件夹: {workflow2_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python integrated_workflow.py <输入目录>")
        sys.exit(1)
    main(os.path.abspath(sys.argv[1]))