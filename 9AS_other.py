#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sort_gff3_by_chr_and_gene.py

功能更新：
1. 自动将文件中所有的 'synthetic' 和 'miniprot' 替换为 'PAP'。
2. 自动删除所有包含 'other' 单词的行。
3. 按染色体和基因起始位置对目录下所有 .gff3 文件重新排序（就地替换）。
"""
import sys
import os
import tempfile
from pathlib import Path
import re

# ---------- 辅助函数 ----------

def extract_seqid_key(seqid: str):
    if seqid is None:
        return (10, seqid or "")
    s = seqid.strip()
    s_tail = s.split("/")[-1]

    # 1) chr + 数字 + 字母 (如 chr01, chr2L)
    m = re.search(r'(?i)chr0*([0-9]+)([A-Za-z]*)$', s_tail)
    if m:
        num = int(m.group(1))
        suffix = m.group(2).upper() if m.group(2) else ""
        return (0, num, suffix, s_tail)

    # 2) chr + 字母 (如 chrX, chrY, chrM)
    m2 = re.search(r'(?i)chr([A-Za-z]+)$', s_tail)
    if m2:
        lab = m2.group(1).upper()
        mapping = {"X": 1000, "Y": 1001, "M": 1002, "MT": 1002}
        mapped = mapping.get(lab, 2000 + (hash(lab) % 1000))
        return (0, mapped, lab, s_tail)

    # 3) 无 chr 前缀 - 提取最后的数字
    m3 = re.search(r'(\d+)(?!.*\d)', s_tail)
    if m3:
        num = int(m3.group(1))
        return (1, num, s_tail)

    # 4) 兜底：纯字符串比较
    return (2, s_tail)

def block_start_coordinate(block_lines):
    for ln in block_lines:
        if ln.startswith("#") or ln.strip() == "":
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) >= 4:
            try:
                return int(cols[3])
            except:
                return float("inf")
    return float("inf")

def block_seqid(block_lines):
    for ln in block_lines:
        if ln.startswith("#") or ln.strip() == "":
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) >= 1:
            return cols[0]
    return None

# ---------- 单个文件处理逻辑 ----------

def process_one_file(path: Path):
    # 读取原始文本
    text = path.read_text(encoding="utf-8")

    # --- 新增步骤 1: 单词替换 ---
    # 将 synthetic 和 miniprot 替换为 PAP
    text = text.replace("synthetic", "PAP").replace("miniprot", "PAP")

    # 分行处理
    raw_lines = text.splitlines(keepends=True)

    # --- 新增步骤 2: 删除含有 'other' 的行 ---
    # 过滤掉任何包含 'other' 字符串的行
    lines = [ln for ln in raw_lines if "other" not in ln]

    if not lines:
        return

    # 1) 提取文件开头的注释 (preamble)
    preamble_lines = []
    idx = 0
    while idx < len(lines) and lines[idx].startswith("#"):
        preamble_lines.append(lines[idx])
        idx += 1

    # 2) 将剩余行解析为 blocks (以 gene 为核心)
    blocks = []
    current_block = None
    i = idx
    while i < len(lines):
        ln = lines[i]
        if ln.strip() == "":
            if current_block is None:
                blocks.append([ln])
            else:
                current_block.append(ln)
            i += 1
            continue
        
        if ln.startswith("#"):
            if current_block is None:
                blocks.append([ln])
            else:
                current_block.append(ln)
            i += 1
            continue

        cols = ln.rstrip("\n").split("\t")
        ftype = cols[2].strip().lower() if len(cols) >= 3 else ""
        
        if ftype == "gene":
            current_block = [ln]
            blocks.append(current_block)
            i += 1
            # 贪婪吸收该 gene 下的所有子结构，直到遇到下一个 gene
            while i < len(lines):
                ln2 = lines[i]
                if ln2.startswith("#") or ln2.strip() == "":
                    current_block.append(ln2)
                    i += 1
                    continue
                cols2 = ln2.rstrip("\n").split("\t")
                f2 = cols2[2].strip().lower() if len(cols2) >= 3 else ""
                if f2 == "gene":
                    break
                current_block.append(ln2)
                i += 1
            current_block = None
        else:
            # 非 gene 类型的独立数据行
            blocks.append([ln])
            i += 1

    # 3) 计算排序键并排序
    block_entries = []
    for bi, block in enumerate(blocks):
        seqid = block_seqid(block)
        start = block_start_coordinate(block)
        chr_key = extract_seqid_key(seqid) if seqid is not None else (99, )
        block_entries.append((chr_key, start, bi, block))

    # 先按染色体(chr_key)，再按起始位置(start)，最后按原索引(bi)稳定排序
    block_entries.sort(key=lambda x: (x[0], x[1], x[2]))

    # 4) 重构并写回文件
    out_lines = preamble_lines + [line for entry in block_entries for line in entry[3]]

    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write("".join(out_lines))
        Path(tmp).replace(path)
    finally:
        if Path(tmp).exists():
            try: Path(tmp).unlink()
            except: pass

def process_dir(dirpath):
    p = Path(dirpath)
    if not p.is_dir():
        print(f"Error: {dirpath} is not a directory", file=sys.stderr)
        return 1
    files = sorted(p.glob("*.gff3"))
    if not files:
        print(f"No .gff3 files under {dirpath}", file=sys.stderr)
        return 1
    for f in files:
        print(f"Processing: {f}")
        try:
            process_one_file(f)
        except Exception as e:
            print(f"Failed {f}: {e}", file=sys.stderr)
    print("Done.")
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 sort_gff3_by_chr_and_gene.py <gff3_directory>", file=sys.stderr)
        sys.exit(2)
    sys.exit(process_dir(sys.argv[1]))
