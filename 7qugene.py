#!/usr/bin/env python3
"""
remove_overlapping_genes_simple.py

在每个染色体（seqid）内按文件顺序遍历 gene：
- 若遇到的 gene 与之前已保留的任意 gene 区间有交集（overlap），则删除该 gene 及其所有后代特征（mRNA/CDS/exon 等）。
- 支持 plain .gff3 和 .gff3.gz。
- 使用 ID= 与 Parent= 字段来追溯层级关系。
"""

import sys
import gzip
from collections import defaultdict

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r", encoding="utf-8")

def parse_attrs(attrstr):
    attrs = {}
    if not attrstr or attrstr == ".":
        return attrs
    parts = attrstr.strip().split(";")
    for p in parts:
        if not p:
            continue
        if "=" in p:
            k, v = p.split("=", 1)
            attrs[k] = v
        elif " " in p:
            k, v = p.split(" ", 1)
            attrs[k] = v
        else:
            attrs[p] = ""
    return attrs

def intervals_overlap(a_start, a_end, b_start, b_end):
    """闭区间重叠判定：a_start <= b_end and b_start <= a_end"""
    if a_start is None or a_end is None or b_start is None or b_end is None:
        return False
    return (a_start <= b_end) and (b_start <= a_end)

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 remove_overlapping_genes_simple.py input.gff3[.gz] > filtered.gff3", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]

    # 读取文件并解析为 entries 列表（保留原始行顺序）
    lines = []
    try:
        fh = open_maybe_gz(path)
    except Exception as e:
        print("Error opening file: {}".format(e), file=sys.stderr)
        sys.exit(2)

    with fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#") or not line.strip():
                lines.append({"raw": line, "is_feature": False})
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                lines.append({"raw": line, "is_feature": False})
                continue
            seqid = cols[0]
            source = cols[1]
            ftype = cols[2]
            try:
                start = int(cols[3])
                end = int(cols[4])
            except:
                start = None
                end = None
            score = cols[5]
            strand = cols[6]
            phase = cols[7]
            attrs_raw = cols[8] if len(cols) > 8 else "."
            attrs = parse_attrs(attrs_raw)
            fid = attrs.get("ID") or attrs.get("id")
            parent = attrs.get("Parent") or attrs.get("parent")
            parents = parent.split(",") if parent else []
            entry = {
                "raw": line,
                "is_feature": True,
                "seqid": seqid,
                "type": ftype,
                "start": start,
                "end": end,
                "attrs_raw": attrs_raw,
                "attrs": attrs,
                "id": fid,
                "parents": parents,
            }
            lines.append(entry)

    # 决定哪些 gene 要删除：若与之前已保留的任一 gene 有交集，则删除该 gene
    removed_gene_ids = set()
    kept_genes_by_seqid = defaultdict(list)  # seqid -> list of (start,end,id) 已保留基因（按文件顺序）

    for e in lines:
        if not e.get("is_feature"):
            continue
        if e["type"].lower() != "gene":
            continue
        seqid = e["seqid"]
        s = e["start"]
        t = e["end"]
        gid = e["id"]
        # 若缺少位置信息或ID，保守地保留
        if s is None or t is None or gid is None:
            kept_genes_by_seqid[seqid].append((s, t, gid))
            continue
        # 检查是否与已保留的任一 gene 有交集
        overlap_any = False
        for ks, kt, kid in kept_genes_by_seqid[seqid]:
            if ks is None or kt is None:
                continue
            if intervals_overlap(s, t, ks, kt):
                overlap_any = True
                break
        if overlap_any:
            removed_gene_ids.add(gid)
        else:
            kept_genes_by_seqid[seqid].append((s, t, gid))

    # 建立 id -> entry index 映射，便于查找 parent 的条目
    id_to_entries = defaultdict(list)
    for idx, e in enumerate(lines):
        if not e.get("is_feature"):
            continue
        if e.get("id"):
            id_to_entries[e["id"]].append(idx)

    # 递归判断某条 entry 是否有被删除的基因作为祖先（Parent 链上是否包含被删除的 gene id）
    removed_cache = {}
    sys.setrecursionlimit(10000)

    def entry_has_removed_ancestor(entry_idx, visited=None):
        if visited is None:
            visited = set()
        if entry_idx in removed_cache:
            return removed_cache[entry_idx]
        if entry_idx in visited:
            removed_cache[entry_idx] = False
            return False
        visited.add(entry_idx)
        e = lines[entry_idx]
        eid = e.get("id")
        if e["type"].lower() == "gene" and eid and eid in removed_gene_ids:
            removed_cache[entry_idx] = True
            return True
        for p in e.get("parents", []):
            if p in removed_gene_ids:
                removed_cache[entry_idx] = True
                return True
            for pi in id_to_entries.get(p, []):
                if entry_has_removed_ancestor(pi, visited):
                    removed_cache[entry_idx] = True
                    return True
        removed_cache[entry_idx] = False
        return False

    # 收集需要删除的 entry 索引（基因本身 + 其子特征）
    entry_idx_to_remove = set()
    for idx, e in enumerate(lines):
        if not e.get("is_feature"):
            continue
        if e["type"].lower() == "gene" and e.get("id") and e["id"] in removed_gene_ids:
            entry_idx_to_remove.add(idx)
            continue
        try:
            if entry_has_removed_ancestor(idx):
                entry_idx_to_remove.add(idx)
        except RecursionError:
            # 若递归过深，保守地不删除该项（避免误删）
            pass

    # 输出保留的行（注释、空行和未被删除的 feature）
    out_lines = []
    for idx, e in enumerate(lines):
        if not e.get("is_feature"):
            out_lines.append(e["raw"])
            continue
        if idx in entry_idx_to_remove:
            continue
        out_lines.append(e["raw"])

    sys.stdout.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
