#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step3_attach_orphan_features_safely.py

- 保持你现有的 orphan 修复 / ID 规则不变（AS_gene... / mRNA / exon / CDS）
- 新增最后一步：对 orphan 区段（marker 之后）进行排序
    1) 按染色体（natural chr 数字优先: chr1, chr2, ...; 非 chrX 按字母序）
    2) 同一染色体内部以 gene 为单位整体移动（gene block 按 gene start 排序）
- 覆盖原文件前会生成 <file>.bak（如不存在）
"""

import sys
import os
import re
from collections import defaultdict

MRNA_RE = re.compile(r'(\d+)\.t\d+')
NUM_RE  = re.compile(r'(\d+)')


def extract_core_num(s):
    """从 attribute 字符串中提取 numeric core（保持原逻辑）"""
    if not s:
        return None
    m = MRNA_RE.search(s)
    if m:
        return m.group(1)
    m = NUM_RE.search(s)
    return m.group(1) if m else None


def ensure_attr_set(attrs, key, value):
    """在 attrs 字符串（a;b=c;...）中设置或替换 key=value，返回新 attrs 字符串"""
    if key + "=" in attrs:
        attrs = re.sub(rf'{key}=[^;]+', f'{key}={value}', attrs)
    else:
        if attrs and attrs.strip():
            attrs = attrs + f";{key}={value}"
        else:
            attrs = f"{key}={value}"
    return attrs


def replace_id_in_attrs(attrs, new_id):
    """把 ID=... 替换为 new_id（若不存在则追加）"""
    if "ID=" in attrs:
        attrs = re.sub(r'ID=[^;]+', f'ID={new_id}', attrs)
    else:
        if attrs and attrs.strip():
            attrs = attrs + f";ID={new_id}"
        else:
            attrs = f"ID={new_id}"
    return attrs


def replace_parent_in_attrs(attrs, new_parent):
    """把 Parent=... 替换为 new_parent（若不存在则追加）"""
    if "Parent=" in attrs:
        attrs = re.sub(r'Parent=[^;]+', f'Parent={new_parent}', attrs)
    else:
        if attrs and attrs.strip():
            attrs = attrs + f";Parent={new_parent}"
        else:
            attrs = f"Parent={new_parent}"
    return attrs


def natural_chr_key(seqid):
    """
    自然排序 key：如果 seqid 中有 chrN (N 为数字) -> (0, N, seqid)
    否则 -> (1, seqid.lower())
    """
    if not seqid:
        return (2, seqid.lower())
    m = re.search(r'chr(\d+)$', seqid, re.IGNORECASE)
    if m:
        return (0, int(m.group(1)), seqid.lower())
    # 有时候 seqid 形如 chr1_extra ——试着从开头匹配
    m2 = re.search(r'chr(\d+)', seqid, re.IGNORECASE)
    if m2:
        return (0, int(m2.group(1)), seqid.lower())
    return (1, seqid.lower())


def process_file(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]

    # 小心替换 marker 字样和 orphan -> AS 字样（保持你之前的替换习惯）
    lines = [
        ln.replace("other_orphan_features", "other_alternative_splicing_features")
          .replace("orphan", "AS")
        for ln in lines
    ]

    # 找到 marker 行
    orphan_idx = None
    for i, ln in enumerate(lines):
        if "other_alternative_splicing_features" in ln:
            orphan_idx = i
            break

    if orphan_idx is None:
        # 没有特殊段，直接返回原样（已替换文本）
        return lines

    head = lines[:orphan_idx + 1]
    tail = lines[orphan_idx + 1:]

    # ---------- 原有修复逻辑（保持不变） ----------
    out_tail = []

    gene_counter = 1
    last_mrna_num = None
    last_mrna_id = None

    exon_counter = defaultdict(int)
    cds_counter = defaultdict(int)

    for line in tail:
        if line.startswith("#") or not line.strip():
            out_tail.append(line)
            continue

        cols = line.split("\t")
        if len(cols) < 9:
            out_tail.append(line)
            continue

        seqid, source, ftype, start, end, score, strand, phase, attrs = cols

        # mRNA: always create synthetic gene and rewrite mRNA ID/Parent
        if ftype == "mRNA":
            gene_id = f"AS_gene{gene_counter}"
            gene_counter += 1
            mrna_id = f"{gene_id}.t01"

            gene_line = "\t".join([
                seqid, "synthetic", "gene",
                start, end, ".", strand, ".",
                f"ID={gene_id}"
            ])
            out_tail.append(gene_line)

            # rewrite mRNA attrs
            attrs = replace_id_in_attrs(attrs, mrna_id)
            attrs = replace_parent_in_attrs(attrs, gene_id)

            cols[8] = attrs
            out_tail.append("\t".join(cols))

            last_mrna_id = mrna_id
            last_mrna_num = extract_core_num(attrs)
            exon_counter[mrna_id] = 0
            cds_counter[mrna_id] = 0
            continue

        # exon / CDS / stop_codon
        if ftype in ("exon", "CDS", "stop_codon"):
            feat_num = extract_core_num(attrs)

            # 若匹配到最近 mRNA (用 numeric core 判定，保持原风格)
            if last_mrna_num and feat_num == last_mrna_num and last_mrna_id:
                # 强制 Parent 指向该 mRNA
                attrs = replace_parent_in_attrs(attrs, last_mrna_id)

                if ftype == "exon":
                    exon_counter[last_mrna_id] += 1
                    new_id = f"{last_mrna_id}-{exon_counter[last_mrna_id]}"
                    attrs = replace_id_in_attrs(attrs, new_id)

                elif ftype == "CDS":
                    cds_counter[last_mrna_id] += 1
                    new_id = f"{last_mrna_id}.cds{cds_counter[last_mrna_id]}"
                    attrs = replace_id_in_attrs(attrs, new_id)

                # stop_codon: 只确保 Parent
                cols[8] = attrs
                out_tail.append("\t".join(cols))
                continue

            # 真正孤立 feature -> 新建 gene + mRNA + feature
            gene_id = f"AS_gene{gene_counter}"
            gene_counter += 1
            mrna_id = f"{gene_id}.t01"

            gene_line = "\t".join([
                seqid, "synthetic", "gene",
                start, end, ".", strand, ".",
                f"ID={gene_id}"
            ])
            mrna_line = "\t".join([
                seqid, "synthetic", "mRNA",
                start, end, ".", strand, ".",
                f"ID={mrna_id};Parent={gene_id}"
            ])

            # rewrite feature Parent -> mrna_id
            attrs = replace_parent_in_attrs(attrs, mrna_id)

            if ftype == "exon":
                exon_counter[mrna_id] = 1
                attrs = replace_id_in_attrs(attrs, f"{mrna_id}-{exon_counter[mrna_id]}")
            elif ftype == "CDS":
                cds_counter[mrna_id] = 1
                attrs = replace_id_in_attrs(attrs, f"{mrna_id}.cds{cds_counter[mrna_id]}")
            # stop_codon: Parent 已替换

            cols[8] = attrs
            out_tail.extend([gene_line, mrna_line, "\t".join(cols)])

            last_mrna_id = None
            last_mrna_num = None
            continue

        # 其它 feature: 原样保留
        out_tail.append(line)

    # ---------- 新增步骤：对 out_tail 进行 chromosome -> gene-block 排序 ----------
    # 将 out_tail 拆分为 prefix_lines（在遇到第一个 gene 前的行，保留顺序）
    # 与 gene_blocks（每个 block 以 gene 行开始，包含随后所有属于该 gene 的行）
    prefix_lines = []
    gene_blocks = []   # list of tuples (seqid, gene_start_int, block_lines)

    cur_block = None
    cur_seqid = None
    cur_start = None
    saw_first_gene = False

    for ln in out_tail:
        if ln.startswith("#") or not ln.strip():
            if not saw_first_gene:
                prefix_lines.append(ln)
            else:
                # comments between blocks — attach to current block if exists, else keep in prefix
                if cur_block is not None:
                    cur_block.append(ln)
                else:
                    prefix_lines.append(ln)
            continue

        cols = ln.split("\t")
        if len(cols) < 9:
            if not saw_first_gene:
                prefix_lines.append(ln)
            else:
                if cur_block is not None:
                    cur_block.append(ln)
                else:
                    prefix_lines.append(ln)
            continue

        ftype = cols[2]
        seqid_line = cols[0]
        if ftype == "gene":
            # start a new block
            if cur_block is not None:
                # push previous
                gene_blocks.append((cur_seqid, cur_start, cur_block))
            saw_first_gene = True
            cur_seqid = seqid_line
            # try parse start
            try:
                cur_start = int(cols[3])
            except:
                cur_start = 10**12
            cur_block = [ln]
        else:
            # non-gene line
            if not saw_first_gene:
                # no gene yet -> put into prefix
                prefix_lines.append(ln)
            else:
                if cur_block is None:
                    # unexpected: attach to prefix
                    prefix_lines.append(ln)
                else:
                    cur_block.append(ln)

    # push last block if any
    if cur_block is not None:
        gene_blocks.append((cur_seqid, cur_start, cur_block))

    # group by seqid
    by_chr = defaultdict(list)
    for seqid_b, start_b, block_b in gene_blocks:
        by_chr[seqid_b].append((start_b, block_b))

    # order seqids naturally
    seqid_list = sorted(by_chr.keys(), key=natural_chr_key)

    sorted_tail = []
    # first keep prefix lines (comments or stray lines before first gene)
    sorted_tail.extend(prefix_lines)

    for seqid_k in seqid_list:
        blocks = by_chr[seqid_k]
        # sort blocks by start
        blocks.sort(key=lambda x: x[0] if isinstance(x[0], int) else 10**12)
        for _, blk in blocks:
            sorted_tail.extend(blk)

    # If there were gene-less lines not captured, append them at end (unlikely)
    # (we already added prefix_lines earlier)
    return head + sorted_tail


def main(gff_dir):
    if not os.path.isdir(gff_dir):
        sys.exit("Not a directory: " + gff_dir)

    for fn in sorted(os.listdir(gff_dir)):
        if not fn.endswith(".gff3"):
            continue
        path = os.path.join(gff_dir, fn)
        out_lines = process_file(path)
        # create backup if not exists
        bak = path + ".bak"
        if not os.path.exists(bak):
            os.rename(path, bak)
        with open(path, "w") as fo:
            for ln in out_lines:
                fo.write(ln + "\n")
        print(f"[OK] {fn} (backup: {os.path.basename(bak)})")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python step3_attach_orphan_features_safely.py <gff_dir>")
    main(sys.argv[1])
