#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
reorder_gff3_exon_cds_inplace.py

- 按 gene block 处理（only before 'other' boundary）
- 删除 stop / stop_codon（仅在 pre-'other' 区域）
- 删除孤立 mRNA（其下没有 exon 和 CDS）【初步删除】
- 全文件级删除：凡是下面没有任何 CDS 的 mRNA（即使有 exon 也认为没有 CDS） -> 删除该 mRNA 及其子层级（此步骤现在对全文执行）
- 若某 gene 下没有任何 mRNA -> 删除该 gene
- 删除后对每个 gene 内剩余 mRNA 顺序紧凑重编号（保留 t 数字宽度），并更新其所有子特征的 Parent/ID
- phase-aware CDS overlap 过滤（同一 mRNA 内）
- overlap collapse（若有重叠，保留先出现的）
- 在输出时为 exon / CDS 重新编号为 <mRNA>.exonN / <mRNA>.CDSN
- 其余非 exon/cds 子特征会把 Parent 更新为新的 mRNA ID（若有）
- 为所有“有 CDS 但在该 mRNA 下没有任何 exon”的 CDS 自动插入一个 exon（坐标完全复制 CDS），插在该 CDS 之前
  （此插入在全文件范围内执行，包含 other 后的部分）
- **保留原始 'other' 边界行**：如果处理链路导致它意外丢失，脚本会把原始边界行按合适位置补回
- 结果就地替换原文件（atomic）
"""

import sys
from pathlib import Path
from collections import OrderedDict
import tempfile
import os
import re

STOP_FEATURES = {"stop", "stop_codon"}

# ---------- utilities ----------

def parse_attrs(s: str):
    d = {}
    for p in s.split(";"):
        if "=" in p:
            k, v = p.split("=", 1)
            d[k] = v
    return d

def format_attrs(d: dict):
    parts = []
    if "ID" in d:
        parts.append(f"ID={d['ID']}")
    if "Parent" in d:
        parts.append(f"Parent={d['Parent']}")
    for k in d:
        if k in ("ID", "Parent"):
            continue
        parts.append(f"{k}={d[k]}")
    return ";".join(parts)

def update_parent_field(parent_val: str, mapping: dict):
    """Parent field may be comma-separated. Replace any parent that matches keys in mapping."""
    parts = [p.strip() for p in parent_val.split(",") if p.strip()]
    newparts = []
    for p in parts:
        newparts.append(mapping.get(p, p))
    return ",".join(newparts)

def overlap_len(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)

def overlap_ratio(a, b):
    ov = overlap_len(a, b)
    ml = min(a[1] - a[0] + 1, b[1] - b[0] + 1)
    return 0.0 if ml == 0 else ov / ml

def parse_phase(cols):
    # cols is a list, phase is column index 7 (0-based)
    try:
        ph = cols[7].strip()
    except Exception:
        return None
    if ph in (".", ""):
        return None
    try:
        p = int(ph)
        return p if p in (0, 1, 2) else None
    except:
        return None

def find_other_boundary(lines):
    # first comment line with text starting 'other' (case-ins), then first data line with 3rd col == 'other'
    for i, ln in enumerate(lines):
        if ln.startswith("#"):
            txt = ln.lstrip("#").strip().lower()
            if txt.startswith("other"):
                return i
    for i, ln in enumerate(lines):
        if ln.startswith("#"):
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) >= 3 and cols[2].strip().lower() == "other":
            return i
    return None

# detect t-number pattern in a mRNA ID, return renumbered id preserving width
_T_RE = re.compile(r'^(.*?t)(\d+)(.*)$', re.IGNORECASE)

def renumber_mrna_id(old_id: str, new_index: int):
    m = _T_RE.match(old_id)
    if m:
        prefix = m.group(1)   # includes 't'
        digits = m.group(2)
        suffix = m.group(3)
        width = len(digits)
        newnum = str(new_index).zfill(width)
        return f"{prefix}{newnum}{suffix}"
    else:
        return f"{old_id}.t{str(new_index).zfill(2)}"

def is_stop_line(ln: str):
    if ln.startswith("#"):
        return False
    cols = ln.rstrip("\n").split("\t")
    if len(cols) < 3:
        return False
    return cols[2].strip().lower() in STOP_FEATURES

# ------------------------- New helper: insert missing exon(s) -------------------------
def insert_exons_for_cds_without_exon(all_lines):
    """
    全文件范围扫描：对于每个 CDS 行，检查该 CDS 的 Parent（可能是逗号分隔多个 parent）。
    对于每个 parent，如果全文件（all_lines）中没有任何 exon 与该 parent 关联，则在 CDS 行之前插入一个 exon 行：
      - exon 的坐标与该 CDS 完全相同
      - exon 的 source 与 CDS 相同
      - exon 的 score '.'，phase '.' （GFF 常规）
      - exon 的 attributes Parent 与 CDS 相同，ID 设为 <first_parent>.autoexonN 临时样式（N 为该 parent 已有 exon 数 + 1）
    返回新的行列表（插入后的）。
    """
    # first scan to find parents that already have exon anywhere
    parent_has_exon = {}
    parent_exon_count = {}  # for temporary numbering
    # Also track existing parent->counts of exons to make unique IDs
    for ln in all_lines:
        if ln.startswith("#"):
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        f = cols[2].strip().lower()
        if f == "exon":
            attrs = parse_attrs(cols[8])
            p = attrs.get("Parent")
            if p:
                for part in [pp.strip() for pp in p.split(",") if pp.strip()]:
                    parent_has_exon[part] = True
                    parent_exon_count[part] = parent_exon_count.get(part, 0) + 1

    new_lines = []
    # Iterate and insert when encountering CDS whose parents lack exons
    for ln in all_lines:
        if ln.startswith("#"):
            new_lines.append(ln)
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 9:
            new_lines.append(ln)
            continue
        f = cols[2].strip().lower()
        if f == "cds":
            attrs = parse_attrs(cols[8])
            pval = attrs.get("Parent")
            if not pval:
                # no Parent: just append
                new_lines.append(ln)
                continue
            parents = [pp.strip() for pp in pval.split(",") if pp.strip()]
            # determine if any parent lacks exon
            need_insert = False
            for pr in parents:
                if not parent_has_exon.get(pr, False):
                    need_insert = True
                    break
            if need_insert:
                # create one exon line preceding this CDS line.
                # We'll set Parent to same comma-separated parents (preserve multi-parent).
                # For ID pick minimal deterministic value: use first parent + .autoexon<N>
                first_parent = parents[0]
                count = parent_exon_count.get(first_parent, 0) + 1
                parent_exon_count[first_parent] = count
                parent_has_exon[first_parent] = True
                # If multiple parents, we won't create separate exon lines; one exon with Parent set to full parent list.
                new_attrs = {"ID": f"{first_parent}.autoexon{count}", "Parent": ",".join(parents)}
                # build exon row: keep seqid, source; type exon; start end same; score '.'; strand same; phase '.'; attrs
                exon_cols = [
                    cols[0],  # seqid
                    cols[1],  # source
                    "exon",
                    cols[3],
                    cols[4],
                    ".",     # score
                    cols[6], # strand
                    ".",     # phase
                    format_attrs(new_attrs)
                ]
                new_lines.append("\t".join(exon_cols) + "\n")
            # then append original CDS line
            new_lines.append(ln)
        else:
            new_lines.append(ln)
    return new_lines

# ------------------------- per-gene initial processing -------------------------

def process_gene_block_initial(block_lines):
    """
    Initial per-gene processing (unchanged logic, but kept here).
    Returns list of lines for this gene block.
    """
    out = []
    n = len(block_lines)
    i = 0

    pre_items = []
    mrnas = []
    current_mrna = None

    while i < n:
        ln = block_lines[i]
        if ln.startswith("#") or ln.strip() == "":
            if current_mrna is None:
                pre_items.append(ln)
            else:
                current_mrna.setdefault("misc_after", []).append(ln)
            i += 1
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 9:
            if current_mrna is None:
                pre_items.append(ln)
            else:
                current_mrna.setdefault("misc_after", []).append(ln)
            i += 1
            continue

        f = cols[2].strip().lower()
        if f == "mrna":
            current_mrna = {
                "mrna_line": ln,
                "mrna_cols": cols,
                "children": [],
                "misc_after": []
            }
            mrnas.append(current_mrna)
            i += 1
            continue

        if current_mrna is not None:
            current_mrna["children"].append((ln, cols))
        else:
            pre_items.append(ln)
        i += 1

    # drop isolated mRNAs (no exon and no cds)
    keep_flags = []
    for mr in mrnas:
        has_exon_or_cds = any(cols[2].strip().lower() in ("exon","cds") for (_, cols) in mr["children"])
        keep_flags.append(has_exon_or_cds)

    # renumber kept mRNAs
    new_ids = {}
    idx = 1
    for keep, mr in zip(keep_flags, mrnas):
        if keep:
            attrs = parse_attrs(mr["mrna_cols"][8])
            oldid = attrs.get("ID") or f"unknown_mrna_{idx}"
            newid = renumber_mrna_id(oldid, idx)
            new_ids[oldid] = newid
            idx += 1

    # helper kept indices
    kept_indices = [i for i, k in enumerate(keep_flags) if k]
    def next_kept(p):
        for idxk in kept_indices:
            if idxk > p:
                return idxk
        return None
    def prev_kept(p):
        for idxk in reversed(kept_indices):
            if idxk < p:
                return idxk
        return None

    # emit pre_items
    for ln in pre_items:
        out.append(ln)

    # process mRNAs
    for i_m, (keep, mr) in enumerate(zip(keep_flags, mrnas)):
        if not keep:
            # reparent other (non-exon/cds) children to neighbor kept mRNA if possible
            for ln, cols in mr["children"]:
                t = cols[2].strip().lower()
                if t in ("exon","cds"):
                    continue
                nk = next_kept(i_m); pk = prev_kept(i_m)
                target_idx = nk if nk is not None else pk
                if target_idx is None:
                    continue
                target_mr = mrnas[target_idx]
                old_target_id = parse_attrs(target_mr["mrna_cols"][8]).get("ID")
                new_target_id = new_ids.get(old_target_id, old_target_id)
                a = parse_attrs(cols[8])
                a["Parent"] = new_target_id
                cols[8] = format_attrs(a)
                target_mr["children"].append((ln, cols))
            continue

        # kept
        attrs = parse_attrs(mr["mrna_cols"][8])
        old_mrna_id = attrs.get("ID")
        new_mrna_id = new_ids.get(old_mrna_id, old_mrna_id)
        attrs["ID"] = new_mrna_id
        mrna_cols = list(mr["mrna_cols"])
        mrna_cols[8] = format_attrs(attrs)
        out.append("\t".join(mrna_cols) + "\n")

        # collect exon/cds by coords and other children
        exon_cds = OrderedDict()
        other_children = []
        for ln, cols in mr["children"]:
            t = cols[2].strip().lower()
            if t in ("exon","cds"):
                try:
                    s = int(cols[3]); e = int(cols[4])
                except:
                    other_children.append((ln, cols))
                    continue
                key = (s,e)
                if key not in exon_cds:
                    exon_cds[key] = {"exon": None, "cds": None}
                if t == "exon":
                    if exon_cds[key]["exon"] is None:
                        exon_cds[key]["exon"] = cols
                else:
                    if exon_cds[key]["cds"] is None:
                        exon_cds[key]["cds"] = cols
            else:
                a = parse_attrs(cols[8])
                a["Parent"] = new_mrna_id
                cols[8] = format_attrs(a)
                other_children.append((ln, cols))

        # build items and filters
        items = []
        for (s,e) in sorted(exon_cds.keys()):
            pair = exon_cds[(s,e)]
            items.append({
                "coord": (s,e),
                "exon": pair["exon"],
                "cds": pair["cds"],
                "phase": parse_phase(pair["cds"]) if pair["cds"] else None
            })

        # phase-aware filter
        kept_items = []
        for iv in items:
            drop = False
            if iv["cds"] is not None:
                for prev in kept_items:
                    if prev["cds"] is not None and overlap_len(prev["coord"], iv["coord"]) > 0:
                        p1 = prev["phase"]; p2 = iv["phase"]
                        if (p1 is not None) and (p2 is not None) and (p1 != p2):
                            drop = True
                            break
            if not drop:
                kept_items.append(iv)

        # collapse overlap
        collapsed = []
        for iv in kept_items:
            if not collapsed:
                collapsed.append(iv)
            else:
                if overlap_ratio(collapsed[-1]["coord"], iv["coord"]) > 0.0:
                    continue
                collapsed.append(iv)

        exon_i = cds_i = 0
        for iv in collapsed:
            if iv["exon"] is not None:
                exon_i += 1
                c = iv["exon"]
                a = parse_attrs(c[8])
                a["ID"] = f"{new_mrna_id}.exon{exon_i}"
                a["Parent"] = new_mrna_id
                c[8] = format_attrs(a)
                out.append("\t".join(c) + "\n")
            if iv["cds"] is not None:
                cds_i += 1
                c = iv["cds"]
                a = parse_attrs(c[8])
                a["ID"] = f"{new_mrna_id}.CDS{cds_i}"
                a["Parent"] = new_mrna_id
                c[8] = format_attrs(a)
                out.append("\t".join(c) + "\n")

        # emit other children
        for ln, cols in other_children:
            out.append("\t".join(cols) + "\n")

    return out

# ------------------------- final global cleanup step (apply to FULL file) -------------------------
def final_global_cleanup(all_lines):
    """
    对 全文 all_lines 做全局清理：
      - 删除任何在全文中没有 CDS 的 mRNA（整个 mRNA block 都删）
      - 删除空 gene（其下没有任何 mRNA）
      - 对保留的每个 gene，按出现顺序紧凑重编号 mRNA（保留 t 宽度），并更新 ID/Parent（包含 other 部分）
    返回清理后的所有行（list of lines）
    """
    n = len(all_lines)
    i = 0
    out = []

    # collect preamble before first gene
    preamble = []
    while i < n:
        ln = all_lines[i]
        if ln.startswith("#"):
            preamble.append(ln); i += 1; continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) >= 3 and cols[2].strip().lower() == "gene":
            break
        preamble.append(ln); i += 1
    out.extend(preamble)

    # iterate gene blocks across whole file
    while i < n:
        if all_lines[i].startswith("#"):
            out.append(all_lines[i]); i += 1; continue
        cols = all_lines[i].rstrip("\n").split("\t")
        if len(cols) < 3 or cols[2].strip().lower() != "gene":
            out.append(all_lines[i]); i += 1; continue

        # collect block
        block_start = i
        i += 1
        while i < n:
            if all_lines[i].startswith("#"):
                i += 1; continue
            c2 = all_lines[i].rstrip("\n").split("\t")
            if len(c2) >= 3 and c2[2].strip().lower() == "gene":
                break
            i += 1
        block_end = i
        block = all_lines[block_start:block_end]

        # parse block into mrnas
        j = 0
        m = len(block)
        gene_line = block[0]
        j = 1
        mrnas = []
        current = None
        pre_items = [gene_line]

        while j < m:
            ln = block[j]
            if ln.startswith("#") or ln.strip() == "":
                if current is None:
                    pre_items.append(ln)
                else:
                    current.setdefault("misc_after", []).append(ln)
                j += 1
                continue
            colsj = ln.rstrip("\n").split("\t")
            if len(colsj) < 9:
                if current is None:
                    pre_items.append(ln)
                else:
                    current.setdefault("misc_after", []).append(ln)
                j += 1
                continue
            f = colsj[2].strip().lower()
            if f == "mrna":
                current = {"mrna_line": ln, "mrna_cols": colsj, "children": [], "misc_after": []}
                mrnas.append(current)
                j += 1
                continue
            if current is not None:
                current["children"].append((ln, colsj))
            else:
                pre_items.append(ln)
            j += 1

        # ----------------- DEDUP LOGIC ADDED HERE -----------------
        # If there are multiple mRNAs with identical mRNA coords (seqid,start,end,strand)
        # then compare their CDS lists:
        #  - if CDS counts differ -> keep all mRNAs
        #  - if CDS counts equal:
        #       - if all lists of CDS coords are identical -> keep only the first mRNA (and drop others)
        #       - else -> keep all mRNAs
        # Dropped duplicates will have their children removed (we assume duplicates are exact and redundant).
        if len(mrnas) > 1:
            # build index -> mRNA coord key
            coord_groups = {}
            for idx_m, mr in enumerate(mrnas):
                try:
                    mr_cols = mr.get("mrna_cols")
                    seqid = mr_cols[0]
                    s = int(mr_cols[3])
                    e = int(mr_cols[4])
                    strand = mr_cols[6]
                    key = (seqid, s, e, strand)
                except Exception:
                    key = None
                coord_groups.setdefault(key, []).append(idx_m)

            duplicate_drop_indices = set()
            for key, idx_list in coord_groups.items():
                if key is None or len(idx_list) <= 1:
                    continue
                # collect CDS coord lists for each mRNA in this group
                cds_lists = []
                counts = []
                for idx_m in idx_list:
                    mr = mrnas[idx_m]
                    cds_coords = []
                    for ln_child, cols_child in mr["children"]:
                        if cols_child[2].strip().lower() == "cds":
                            try:
                                cs = int(cols_child[3]); ce = int(cols_child[4])
                                cds_coords.append((cs, ce))
                            except:
                                pass
                    cds_coords_sorted = sorted(cds_coords)
                    cds_lists.append(cds_coords_sorted)
                    counts.append(len(cds_coords_sorted))
                # if counts differ -> skip (do nothing)
                if len(set(counts)) != 1:
                    continue
                # all counts same; check if all cds_lists identical
                all_same = True
                first = cds_lists[0]
                for other in cds_lists[1:]:
                    if other != first:
                        all_same = False
                        break
                if all_same:
                    # drop all but first index in idx_list
                    # keep the earliest appearing mRNA (smallest idx_m)
                    keep_idx = min(idx_list)
                    for idx_m in idx_list:
                        if idx_m != keep_idx:
                            duplicate_drop_indices.add(idx_m)
            # mark dropped duplicates by emptying their children and tagging
            for di in duplicate_drop_indices:
                mrnas[di]["_dropped_duplicate"] = True
                # remove children to avoid reparenting or emitting them later
                mrnas[di]["children"] = []
        # ----------------- END DEDUP LOGIC -----------------

        # Decide deletion: any mRNA that has NO CDS among its children (in the entire file) should be deleted.
        keep_flags = []
        for mr in mrnas:
            # if this mRNA was marked as dropped duplicate, treat it as having no CDS (so it will be removed)
            if mr.get("_dropped_duplicate"):
                keep_flags.append(False)
                continue
            has_cds = any(cols[2].strip().lower() == "cds" for (_, cols) in mr["children"])
            keep_flags.append(has_cds)

        # if no mRNA kept -> delete whole gene block
        if not any(keep_flags):
            continue

        # renumber kept mRNAs
        new_ids = {}
        idx = 1
        for keep, mr in zip(keep_flags, mrnas):
            if keep:
                attrs = parse_attrs(mr["mrna_cols"][8])
                oldid = attrs.get("ID") or f"unknown_mrna_{idx}"
                newid = renumber_mrna_id(oldid, idx)
                new_ids[oldid] = newid
                idx += 1

        # reparent non-exon/cds children of dropped mRNAs to nearest kept if possible
        kept_indices = [ii for ii, k in enumerate(keep_flags) if k]
        def next_kept(pos):
            for kk in kept_indices:
                if kk > pos:
                    return kk
            return None
        def prev_kept(pos):
            for kk in reversed(kept_indices):
                if kk < pos:
                    return kk
            return None

        for pos, (keep, mr) in enumerate(zip(keep_flags, mrnas)):
            if not keep:
                # If this mRNA was a duplicate we dropped intentionally, do NOT reparent its children:
                # we already removed its children in dedup step to avoid duplication.
                if mr.get("_dropped_duplicate"):
                    continue
                for ln, cols in mr["children"]:
                    t = cols[2].strip().lower()
                    if t in ("exon","cds"):
                        continue
                    nk = next_kept(pos); pk = prev_kept(pos)
                    target = nk if nk is not None else pk
                    if target is None:
                        continue
                    target_mr = mrnas[target]
                    old_target_id = parse_attrs(target_mr["mrna_cols"][8]).get("ID")
                    new_target_id = new_ids.get(old_target_id, old_target_id)
                    a = parse_attrs(cols[8])
                    a["Parent"] = new_target_id
                    cols[8] = format_attrs(a)
                    target_mr["children"].append((ln, cols))

        # emit gene_line
        out.append(gene_line)

        # emit each kept mrna, renumber and update children
        for pos, (keep, mr) in enumerate(zip(keep_flags, mrnas)):
            if not keep:
                continue
            old_mrna_id = parse_attrs(mr["mrna_cols"][8]).get("ID")
            new_mrna_id = new_ids.get(old_mrna_id, old_mrna_id)

            # emit updated mRNA line
            attrs = parse_attrs(mr["mrna_cols"][8])
            attrs["ID"] = new_mrna_id
            new_mrna_cols = list(mr["mrna_cols"])
            new_mrna_cols[8] = format_attrs(attrs)
            out.append("\t".join(new_mrna_cols) + "\n")

            # collect exon/cds and other children, update their Parent/ID
            exon_cds = OrderedDict()
            other_children = []
            for ln, cols in mr["children"]:
                t = cols[2].strip().lower()
                if t in ("exon","cds"):
                    try:
                        s = int(cols[3]); e = int(cols[4])
                    except:
                        other_children.append((ln, cols)); continue
                    key = (s,e)
                    if key not in exon_cds:
                        exon_cds[key] = {"exon": None, "cds": None}
                    if t == "exon":
                        if exon_cds[key]["exon"] is None:
                            exon_cds[key]["exon"] = cols
                    else:
                        if exon_cds[key]["cds"] is None:
                            exon_cds[key]["cds"] = cols
                else:
                    a = parse_attrs(cols[8])
                    pval = a.get("Parent")
                    if pval:
                        newp = update_parent_field(pval, {old_mrna_id: new_mrna_id})
                        a["Parent"] = newp
                    else:
                        a["Parent"] = new_mrna_id
                    cols[8] = format_attrs(a)
                    other_children.append((ln, cols))

            # build items
            items = []
            for (s,e) in sorted(exon_cds.keys()):
                pair = exon_cds[(s,e)]
                items.append({
                    "coord": (s,e),
                    "exon": pair["exon"],
                    "cds": pair["cds"],
                    "phase": parse_phase(pair["cds"]) if pair["cds"] else None
                })

            # phase-aware filter
            kept_items = []
            for iv in items:
                drop = False
                if iv["cds"] is not None:
                    for prev in kept_items:
                        if prev["cds"] is not None and overlap_len(prev["coord"], iv["coord"]) > 0:
                            p1 = prev["phase"]; p2 = iv["phase"]
                            if (p1 is not None) and (p2 is not None) and (p1 != p2):
                                drop = True
                                break
                if not drop:
                    kept_items.append(iv)

            # collapse overlap
            collapsed = []
            for iv in kept_items:
                if not collapsed:
                    collapsed.append(iv)
                else:
                    if overlap_ratio(collapsed[-1]["coord"], iv["coord"]) > 0.0:
                        continue
                    collapsed.append(iv)

            # output interleaved exon/CDS renumbered
            exon_i = cds_i = 0
            for iv in collapsed:
                if iv["exon"] is not None:
                    exon_i += 1
                    c = iv["exon"]
                    a = parse_attrs(c[8])
                    a["ID"] = f"{new_mrna_id}.exon{exon_i}"
                    a["Parent"] = new_mrna_id
                    c[8] = format_attrs(a)
                    out.append("\t".join(c) + "\n")
                if iv["cds"] is not None:
                    cds_i += 1
                    c = iv["cds"]
                    a = parse_attrs(c[8])
                    a["ID"] = f"{new_mrna_id}.CDS{cds_i}"
                    a["Parent"] = new_mrna_id
                    c[8] = format_attrs(a)
                    out.append("\t".join(c) + "\n")

            # emit other children
            for ln, cols in other_children:
                out.append("\t".join(cols) + "\n")

    return out

# ------------------------- top-level processing -------------------------

def process_file_inplace(path_in: str):
    p = Path(path_in)
    text = p.read_text()
    lines = text.splitlines(keepends=True)

    # record original other-boundary line(s) (first matching). We'll preserve this exact string if possible.
    bound = find_other_boundary(lines)
    original_boundary_line = None
    original_post_first = None
    if bound is not None:
        original_boundary_line = lines[bound]
        # remember next post line (if any) to try to insert boundary back before it if missing later
        if bound + 1 < len(lines):
            original_post_first = lines[bound + 1]

    # split pre/post around boundary (post includes boundary line)
    if bound is None:
        pre = lines[:]
        post = []
    else:
        pre = lines[:bound]
        post = lines[bound:]  # includes boundary line

    # drop stop lines in pre (preserve comments and other lines)
    filtered_pre = [ln for ln in pre if not is_stop_line(ln)]
    pre = filtered_pre

    # split pre into gene blocks, process each block initial pass
    i = 0
    n = len(pre)
    output_pre = []

    # collect preamble before first gene
    while i < n:
        ln = pre[i]
        if ln.startswith("#"):
            output_pre.append(ln); i += 1; continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) >= 3 and cols[2].strip().lower() == "gene":
            break
        output_pre.append(ln); i += 1

    # iterate gene blocks
    while i < n:
        if pre[i].startswith("#"):
            output_pre.append(pre[i]); i += 1; continue
        cols = pre[i].rstrip("\n").split("\t")
        if len(cols) >= 3 and cols[2].strip().lower() == "gene":
            start = i; i += 1
            while i < n:
                if pre[i].startswith("#"):
                    i += 1; continue
                c2 = pre[i].rstrip("\n").split("\t")
                if len(c2) >= 3 and c2[2].strip().lower() == "gene":
                    break
                i += 1
            end = i
            block = pre[start:end]
            processed_block = process_gene_block_initial(block)
            output_pre.extend(processed_block)
            continue
        else:
            output_pre.append(pre[i]); i += 1

    # **关键点**：合并 pre（初步处理结果）和 post（原始其他部分）
    combined = output_pre + post

    # 新增一步：在全文件范围内为“有 CDS 但该 mRNA 没有任何 exon”的 CDS 插入 exon（插在该 CDS 之前）
    combined_with_exons = insert_exons_for_cds_without_exon(combined)

    # 再做 final 全局清理（该函数会删除无 CDS 的 mRNA、删除空 gene，并重编号等）
    cleaned_all = final_global_cleanup(combined_with_exons)

    # final result
    final_lines = cleaned_all

    # Ensure original boundary line is present — if missing, try to insert it back before original_post_first (if found),
    # otherwise append it near the split (best-effort).
    if original_boundary_line is not None and original_boundary_line not in final_lines:
        inserted = False
        if original_post_first is not None:
            # try find original_post_first in final_lines
            try:
                idx = final_lines.index(original_post_first)
                final_lines.insert(idx, original_boundary_line)
                inserted = True
            except ValueError:
                inserted = False
        if not inserted:
            # fallback: append boundary near middle — we append at the end of the pre-region:
            # find first gene line; insert boundary before first gene if exists, else append at top
            placed = False
            for k, ln in enumerate(final_lines):
                if not ln.startswith("#"):
                    cols = ln.rstrip("\n").split("\t")
                    if len(cols) >= 3 and cols[2].strip().lower() == "gene":
                        final_lines.insert(k, original_boundary_line)
                        placed = True
                        break
            if not placed:
                # last resort append at end
                final_lines.append(original_boundary_line)

    # atomic write
    dirn = p.parent
    fd, tmp = tempfile.mkstemp(dir=str(dirn), prefix=p.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write("".join(final_lines))
        Path(tmp).replace(p)
    finally:
        if Path(tmp).exists():
            try:
                Path(tmp).unlink()
            except:
                pass

def process_path(path_str: str):
    p = Path(path_str)
    if p.is_dir():
        files = sorted(p.glob("*.gff3"))
        if not files:
            print(f"No .gff3 files under {p}", file=sys.stderr)
            return
        for f in files:
            print(f"Processing (in-place): {f}")
            process_file_inplace(str(f))
        print("Done.")
    else:
        p = Path(path_str)
        if not p.exists() or p.suffix.lower() != ".gff3":
            print("Input must be a directory or a .gff3 file", file=sys.stderr)
            sys.exit(1)
        print(f"Processing (in-place): {p}")
        process_file_inplace(str(p))
        print("Done.")

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 reorder_gff3_exon_cds_inplace.py <gff3 file or directory>", file=sys.stderr)
        sys.exit(1)
    process_path(sys.argv[1])

if __name__ == "__main__":
    main()
