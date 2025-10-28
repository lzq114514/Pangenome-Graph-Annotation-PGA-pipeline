#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_miniprot_and_group.py

合并两个脚本：
 - 原 batch_miniprot_and_strictfilter.py 的 sbatch 提交逻辑（每个 .fa 生成一个作业）
 - fast_group_and_renumber_rename_feats.py 的分组/重命名输出逻辑（保证 gene 在其包含的特征前面，并重写非-gene/mRNA 特征 ID/Parent）

用法举例:
    python3 batch_miniprot_and_group.py --input-dir /path/to/genomes_dir --protein /path/to/miniprotzhushi.fa --threads 16 --group-workers 2

注意:
 - 需要系统上可用 miniprot 和 sbatch（如果你希望提交 sbatch）。
 - 该脚本会在每个样本目录下写入三个 helper 脚本并提交 sbatch 作业。
"""
from __future__ import annotations
import os
import sys
import glob
import subprocess
from pathlib import Path
from typing import List

# ---------------- utilities ----------------
def find_genome_files(input_dir: str) -> List[str]:
    exts = ("*.fna", "*.fa", "*.fasta", "*.FNA", "*.FA", "*.FASTA")
    files = []
    for e in exts:
        files.extend(sorted(glob.glob(os.path.join(input_dir, e))))
    return files

def safe_mkdir(path: str):
    os.makedirs(path, exist_ok=True)

# ---------------- SBATCH template ----------------
SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --error=%j.err
#SBATCH --output=%j.out

set -euo pipefail
cd "{job_dir}"

echo "Job start: $(date)"
echo "Genome: {genome}"
echo "Protein: {protein}"
echo "Threads: {threads}"
echo "Index file: {index_file}"
echo "Raw GFF: {raw_gff}"
echo "Fixed step2 GFF: {fixed_step2}"
echo "Strict GFF: {strict_gff}"
echo "Final GFF3: {final_gff}"

# Build index
miniprot -t{threads} -d "{index_file}" "{genome}"

# Align and write raw GFF
miniprot -t{threads} --gff "{index_file}" "{protein}" > "{raw_gff}"

# 2) normalize + ensure each mRNA has a gene parent
python3 "{fix_script}" "{raw_gff}" "{fixed_step2}"

# 3) strict filtering (identical/overlap/delete, containment->reparent)
python3 "{strict_script}" "{fixed_step2}" "{strict_gff}"

# 4) group, renumber and rename features; ensure gene before contained features
python3 "{group_script}" "{strict_gff}" -o "{final_gff}" --tmpdir "{job_dir}/tmp" --workers {group_workers}

echo "Job done: $(date)"
"""

# ---------------- helper scripts contents ----------------
# fix_and_ensure.py (from your first script)
FIX_AND_ENSURE_PY = r'''#!/usr/bin/env python3
# fix_and_ensure.py
# 读取 raw miniprot gff（按 tab），生成规范化 GFF，并保证每个 mRNA 有 gene 父本。
import sys, re

def parse_line_keep_attr(line):
    # 按 tab 分割，保持第9列（attrs）原样（即使其中有空格）
    cols = line.rstrip("\n").split("\t")
    # 如果列数<9，补齐，attrs为"."
    while len(cols) < 9:
        cols.append('.')
    return cols

def fix_gff(in_path, out_path):
    gene_count = 1
    exon_count = 1
    current_parent = None
    with open(in_path, 'r', encoding='utf-8', errors='replace') as fin, open(out_path, 'w', encoding='utf-8') as fout:
        fout.write("##gff-version 3\n")
        for raw in fin:
            if not raw.strip() or raw.startswith("#") or raw.startswith("##PAF"):
                continue
            cols = parse_line_keep_attr(raw)
            seqid = cols[0]
            source = cols[1] if len(cols) > 1 else '.'
            feature = cols[2] if len(cols) > 2 else '.'
            start = cols[3] if len(cols) > 3 else '.'
            end = cols[4] if len(cols) > 4 else '.'
            score = cols[5] if len(cols) > 5 else '.'
            strand = cols[6] if len(cols) > 6 else '.'
            phase = cols[7] if len(cols) > 7 else '.'
            attr = cols[8] if len(cols) > 8 else '.'

            # 规范化 mRNA -> 先写 gene，再写 mRNA，确保 Parent 指向 gene
            if feature == "mRNA":
                exon_count = 1
                m = re.search(r'ID=([^;\s]+)', attr)
                if m:
                    mrna_id = m.group(1)
                else:
                    mrna_id = "MP" + str(gene_count).zfill(6)
                    # 保持 attr 原样，但加 ID=...
                    if attr == ".":
                        attr = f"ID={mrna_id}"
                    else:
                        attr = f"ID={mrna_id};" + attr
                gene_id = "gene-" + mrna_id
                gene_attr = f"ID={gene_id};Name={gene_id}"
                fout.write("\t".join([seqid, source, "gene", start, end, score, strand, ".", gene_attr]) + "\n")
                # 确保 mRNA 的 Parent 指向 gene_id
                if 'Parent=' in attr:
                    attr = re.sub(r'Parent=[^;\s]+', "Parent=" + gene_id, attr)
                else:
                    if 'ID=' in attr:
                        attr = re.sub(r'ID=[^;\s]+', f"ID={mrna_id};Parent={gene_id}", attr)
                    else:
                        attr = f"ID={mrna_id};Parent={gene_id};" + (attr if attr != '.' else "")
                fout.write("\t".join([seqid, source, "mRNA", start, end, score, strand, ".", attr]) + "\n")
                current_parent = mrna_id
                gene_count += 1
                continue

            # exon under current_parent
            if feature == "exon" and current_parent:
                if 'ID=' not in attr or attr == ".":
                    exon_id = f"exon-{current_parent}-{exon_count}"
                    exon_count += 1
                    attr = f"ID={exon_id};Parent={current_parent}"
                else:
                    if 'Parent=' in attr:
                        attr = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                    else:
                        attr = attr + f";Parent={current_parent}"
                fout.write("\t".join([seqid, source, "exon", start, end, score, strand, phase, attr]) + "\n")
                continue

            # CDS -> generate exon then CDS (CDS keep Parent)
            if feature == "CDS" and current_parent:
                exon_id = f"exon-{current_parent}-{exon_count}"
                exon_count += 1
                exon_attr = f"ID={exon_id};Parent={current_parent}"
                fout.write("\t".join([seqid, source, "exon", start, end, ".", strand, ".", exon_attr]) + "\n")
                if 'Parent=' in attr:
                    attr_cds = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                else:
                    if attr == ".":
                        attr_cds = f"Parent={current_parent}"
                    else:
                        attr_cds = f"Parent={current_parent};" + attr
                fout.write("\t".join([seqid, source, "CDS", start, end, score, strand, phase, attr_cds]) + "\n")
                continue

            # stop_codon
            if feature == "stop_codon" and current_parent:
                if 'Parent=' in attr:
                    attr_stop = re.sub(r'Parent=[^;\s]+', f"Parent={current_parent}", attr)
                else:
                    if attr == ".":
                        attr_stop = f"Parent={current_parent}"
                    else:
                        attr_stop = f"Parent={current_parent};" + attr
                fout.write("\t".join([seqid, source, "stop_codon", start, end, score, strand, phase, attr_stop]) + "\n")
                continue

            # 其它直接输出（保持第9列原样，不做空格拆分）
            fout.write("\t".join([seqid, source, feature, start, end, score, strand, phase, attr]) + "\n")

def ensure_genes(in_path, out_path):
    # 如果某些 mRNA 没有 Parent（gene），在其之前插入 gene 行并修改 Parent 指向。
    with open(in_path, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()

    out_lines = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("#"):
            out_lines.append(ln)
            i += 1
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 9:
            out_lines.append(ln)
            i += 1
            continue
        ftype = cols[2].lower()
        attr = cols[8]
        if ftype == "mrna":
            parent_match = re.search(r'Parent=([^;\s]+)', attr)
            id_match = re.search(r'ID=([^;\s]+)', attr)
            parent = parent_match.group(1) if parent_match else None
            mid = id_match.group(1) if id_match else None
            if not parent:
                if mid:
                    new_gene_id = "gene-" + mid
                else:
                    new_gene_id = "gene_unspecified_" + str(i)
                gene_attr = "ID=" + new_gene_id + ";Name=" + new_gene_id
                gene_line = "\t".join([cols[0], cols[1], "gene", cols[3], cols[4], cols[5], cols[6], cols[7], gene_attr]) + "\n"
                out_lines.append(gene_line)
                # update mRNA's Parent
                if 'Parent=' in attr:
                    attr = re.sub(r'Parent=[^;\s]+', "Parent=" + new_gene_id, attr)
                else:
                    if 'ID=' in attr:
                        attr = re.sub(r'ID=[^;\s]+', "ID=" + (mid if mid else "unknown") + ";Parent=" + new_gene_id, attr)
                    else:
                        attr = "ID=" + (mid if mid else "unknown") + ";Parent=" + new_gene_id + ";" + (attr if attr != '.' else "")
                cols[8] = attr
                out_lines.append("\t".join(cols) + "\n")
                i += 1
                continue
        out_lines.append(ln)
        i += 1

    with open(out_path, "w", encoding="utf-8") as outfh:
        outfh.writelines(out_lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 fix_and_ensure.py raw.gff fixed.step2.gff", file=sys.stderr)
        sys.exit(1)
    raw = sys.argv[1]
    step2 = sys.argv[2]
    fix_gff(raw, step2 + ".tmp")
    ensure_genes(step2 + ".tmp", step2)
    # remove tmp
    import os
    try:
        os.remove(step2 + ".tmp")
    except:
        pass
'''

# ---------------- strict_filter_merge.py (from your first script) ----------------
STRICT_FILTER_MERGE_PY = r'''#!/usr/bin/env python3
# strict_filter_merge.py
# (content identical to what you provided)
import sys, gzip, re
from collections import defaultdict, deque

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r", encoding="utf-8")

def parse_attrs(attrstr):
    d = {}
    if not attrstr or attrstr == ".":
        return d
    for part in attrstr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k,v = part.split("=",1)
            d[k] = v
        else:
            d[part] = ""
    return d

def build_attrstr(attrs):
    if not attrs:
        return "."
    return ";".join([f"{k}={v}" if v!="" else k for k,v in attrs.items()])

def intervals_overlap(a_start,a_end,b_start,b_end):
    if None in (a_start,a_end,b_start,b_end):
        return False
    return (a_start <= b_end) and (b_start <= a_end)

def is_strictly_contained(child_s, child_e, parent_s, parent_e):
    if None in (child_s, child_e, parent_s, parent_e):
        return False
    return (child_s > parent_s) and (child_e < parent_e)

def strip_t_suffix(mid):
    return re.sub(r'(\.t\d+)$', '', mid)

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 strict_filter_merge.py input.gff out.gff", file=sys.stderr)
        sys.exit(1)
    in_path = sys.argv[1]
    out_path = sys.argv[2]

    entries = []
    try:
        fh = open_maybe_gz(in_path)
    except Exception as e:
        print("Error opening file:", e, file=sys.stderr)
        sys.exit(2)

    with fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#") or not line.strip():
                entries.append({"is_feature": False, "raw": line})
                continue
            cols = line.split("\t")
            while len(cols) < 9:
                cols.append(".")
            seqid = cols[0]
            ftype = cols[2]
            try:
                start = int(cols[3])
                end = int(cols[4])
            except:
                start = None
                end = None
            attrstr = cols[8]
            attrs = parse_attrs(attrstr)
            fid = attrs.get("ID") or attrs.get("id")
            parent_raw = attrs.get("Parent") or attrs.get("parent") or ""
            parents = parent_raw.split(",") if parent_raw else []
            entries.append({
                "is_feature": True,
                "raw": line,
                "cols": cols,
                "seqid": seqid,
                "type": ftype,
                "start": start,
                "end": end,
                "attrs": attrs,
                "id": fid,
                "parents": parents,
            })

    # id -> indices
    id_to_idxs = defaultdict(list)
    parent_to_children = defaultdict(list)
    for i,e in enumerate(entries):
        if not e.get("is_feature"):
            continue
        if e.get("id"):
            id_to_idxs[e["id"]].append(i)
        for p in e.get("parents", []):
            parent_to_children[p].append(i)

    genes_by_seq = defaultdict(list)
    for i,e in enumerate(entries):
        if not e.get("is_feature"): continue
        if e["type"].lower() == "gene":
            genes_by_seq[e["seqid"]].append((i, e.get("id"), e.get("start"), e.get("end")))

    kept_genes_by_seq = defaultdict(list)
    removed_gene_ids = set()
    contained_map = {}

    for seq, glist in genes_by_seq.items():
        for idx,gid,s,e in glist:
            if gid is None or s is None or e is None:
                kept_genes_by_seq[seq].append((idx,gid,s,e))
                continue
            identical = False
            contained_by = None
            overlaps_any = False
            overlap_with_kid = None
            for (kidx, kid, ks, kt) in kept_genes_by_seq[seq]:
                if ks is None or kt is None:
                    continue
                if s == ks and e == kt:
                    identical = True
                    break
                if is_strictly_contained(s,e,ks,kt):
                    contained_by = kid
                    break
                if intervals_overlap(s,e,ks,kt):
                    overlaps_any = True
                    overlap_with_kid = kid
            if identical:
                removed_gene_ids.add(gid)
                continue
            if contained_by:
                contained_map[gid] = contained_by
                continue
            if overlaps_any:
                contained_map[gid] = overlap_with_kid
                continue
            kept_genes_by_seq[seq].append((idx,gid,s,e))

    need_synth_mrna_for = {}
    synth_counter = 0
    def make_synth_id(parent_gene_id):
        nonlocal synth_counter
        synth_counter += 1
        return f"{parent_gene_id}.auto_mrna{synth_counter}"

    for child_gid, parent_gid in contained_map.items():
        child_idxs = parent_to_children.get(child_gid, [])
        has_mrna = False
        has_other = False
        for ci in child_idxs:
            ce = entries[ci]
            if ce["type"].lower() == "mrna":
                has_mrna = True
            else:
                if ce["type"].lower() == "gene":
                    continue
                has_other = True
        if not has_mrna and has_other:
            need_synth_mrna_for[child_gid] = None

    to_remove_entry_idxs = set()
    from collections import deque
    q = deque()
    for gid in removed_gene_ids:
        for gi in id_to_idxs.get(gid, []):
            ge = entries[gi]
            if ge["type"].lower() == "gene":
                q.append(gi)
                to_remove_entry_idxs.add(gi)
    while q:
        cur = q.popleft()
        cur_e = entries[cur]
        cur_id = cur_e.get("id")
        if not cur_id:
            continue
        for child_idx in parent_to_children.get(cur_id, []):
            if child_idx not in to_remove_entry_idxs:
                to_remove_entry_idxs.add(child_idx)
                q.append(child_idx)

    for child_gid in contained_map.keys():
        for gi in id_to_idxs.get(child_gid, []):
            ge = entries[gi]
            if ge["type"].lower() == "gene":
                to_remove_entry_idxs.add(gi)

    gene_to_mrna_idxs = defaultdict(list)
    for gid, child_idxs in parent_to_children.items():
        for ci in child_idxs:
            ce = entries[ci]
            if ce["type"].lower() == "mrna":
                gene_to_mrna_idxs[gid].append(ci)

    for child_gid, parent_gid in contained_map.items():
        if child_gid == parent_gid:
            continue
        child_mrnas = gene_to_mrna_idxs.get(child_gid, [])
        if child_mrnas:
            gene_to_mrna_idxs[parent_gid].extend(child_mrnas)

    gene_info = {}
    for seq, glist in genes_by_seq.items():
        for idx, gid, s, e in glist:
            if gid:
                gene_info[gid] = (idx, s, e)

    parent_new_bounds = {}
    for gid, mrna_idx_list in gene_to_mrna_idxs.items():
        if gid not in gene_info:
            continue
        gidx, gs, ge = gene_info[gid]
        min_s = None
        max_e = None
        for mi in mrna_idx_list:
            m_ent = entries[mi]
            ms = m_ent.get("start")
            me = m_ent.get("end")
            if ms is None or me is None:
                continue
            if min_s is None or (ms < min_s):
                min_s = ms
            if max_e is None or (me > max_e):
                max_e = me
        if min_s is None:
            min_s = gs
        if max_e is None:
            max_e = ge
        new_s = min(gs if gs is not None else min_s, min_s)
        new_e = max(ge if ge is not None else max_e, max_e)
        parent_new_bounds[gid] = (new_s, new_e)

    old_mrna_to_new = {}
    gene_to_mrna_neworder = {}
    for gid, mrna_idxs in gene_to_mrna_idxs.items():
        unique_idxs = sorted(set(mrna_idxs), key=lambda ii: (entries[ii].get("start") if entries[ii].get("start") is not None else 10**18, ii))
        if not unique_idxs:
            continue
        first_idx = unique_idxs[0]
        first_ent = entries[first_idx]
        first_old_id = first_ent.get("id") or f"{gid}.r1"
        base = strip_t_suffix(first_old_id)
        enumerated = []
        for i, mi in enumerate(unique_idxs, start=1):
            new_id = f"{base}.t{int(i):02d}"
            old_id = entries[mi].get("id") or f"{base}.old{i}"
            old_mrna_to_new[old_id] = new_id
            enumerated.append((mi, old_id, new_id))
        gene_to_mrna_neworder[gid] = enumerated

    for gid, (gidx, gs, ge) in gene_info.items():
        if gid in gene_to_mrna_neworder:
            continue
        mrna_idxs = gene_to_mrna_idxs.get(gid, [])
        if not mrna_idxs:
            continue
        unique_idxs = sorted(set(mrna_idxs), key=lambda ii: (entries[ii].get("start") if entries[ii].get("start") is not None else 10**18, ii))
        if not unique_idxs:
            continue
        first_idx = unique_idxs[0]
        first_old_id = entries[first_idx].get("id") or f"{gid}.r1"
        base = strip_t_suffix(first_old_id)
        enumerated = []
        for i, mi in enumerate(unique_idxs, start=1):
            new_id = f"{base}.t{int(i):02d}"
            old_id = entries[mi].get("id") or f"{base}.old{i}"
            old_mrna_to_new[old_id] = new_id
            enumerated.append((mi, old_id, new_id))
        gene_to_mrna_neworder[gid] = enumerated

    mrna_old_to_parentgene = {}
    for gid, enumerated in gene_to_mrna_neworder.items():
        for (mi, old_id, new_id) in enumerated:
            mrna_old_to_parentgene[old_id] = gid
            mrna_old_to_parentgene[new_id] = gid

    out_entries = []
    emitted_synth = {}
    parent_map = contained_map

    for idx,e in enumerate(entries):
        if not e.get("is_feature"):
            out_entries.append(e["raw"])
            continue
        if idx in to_remove_entry_idxs:
            continue

        cols = list(e["cols"])
        etype = e["type"].lower()
        attrs = dict(e["attrs"])

        if etype == "gene":
            gid = e.get("id")
            if gid in parent_new_bounds:
                new_s, new_e = parent_new_bounds[gid]
                cols[3] = str(new_s)
                cols[4] = str(new_e)
                cols[8] = build_attrstr(attrs)
                out_entries.append("\t".join(cols))
                continue
            else:
                out_entries.append(e["raw"])
                continue

        if etype == "mrna":
            old_id = e.get("id")
            orig_parents = e.get("parents") or []
            assigned_parent = None
            for p in orig_parents:
                if p in parent_map:
                    assigned_parent = parent_map[p]
                    break
                else:
                    assigned_parent = p
                    break
            new_id = None
            if old_id in old_mrna_to_new:
                new_id = old_mrna_to_new[old_id]
            else:
                for gid, enumerated in gene_to_mrna_neworder.items():
                    for (mi, oldi, newi) in enumerated:
                        if mi == idx:
                            new_id = newi
                            old_id = oldi
                            assigned_parent = gid
                            break
                    if new_id:
                        break
            if new_id is None:
                if 'Parent' in attrs:
                    pval = attrs['Parent']
                    p_list = pval.split(",")
                    new_p_list = []
                    for p in p_list:
                        if p in parent_map:
                            new_p_list.append(parent_map[p])
                        else:
                            new_p_list.append(p)
                    attrs['Parent'] = ",".join(new_p_list)
                cols[8] = build_attrstr(attrs)
                out_entries.append("\t".join(cols))
                continue

            attrs['ID'] = new_id
            if assigned_parent:
                attrs['Parent'] = assigned_parent
            cols[8] = build_attrstr(attrs)
            out_entries.append("\t".join(cols))
            continue

        if 'Parent' in attrs or 'parent' in attrs:
            pkey = 'Parent' if 'Parent' in attrs else 'parent'
            pval = attrs.get(pkey, "")
            p_list = pval.split(",") if pval else []
            new_parents = []
            changed = False
            for p in p_list:
                if p in parent_map:
                    par = parent_map[p]
                    enumerated = gene_to_mrna_neworder.get(par, [])
                    if enumerated:
                        new_par_id = enumerated[0][2]
                        new_parents.append(new_par_id)
                        changed = True
                    else:
                        new_parents.append(par)
                        changed = True
                elif p in old_mrna_to_new:
                    new_parents.append(old_mrna_to_new[p])
                    changed = True
                else:
                    new_parents.append(p)
            if changed:
                attrs[pkey] = ",".join(new_parents)
                cols[8] = build_attrstr(attrs)
                out_entries.append("\t".join(cols))
                continue
            else:
                out_entries.append(e["raw"])
                continue
        else:
            out_entries.append(e["raw"])

    with open(out_path, "w", encoding="utf-8") as outf:
        wrote_version = False
        for ln in out_entries:
            if ln.strip().lower().startswith("##gff-version") and not wrote_version:
                outf.write(ln.rstrip("\n") + "\n")
                wrote_version = True
        if not wrote_version:
            outf.write("##gff-version 3\n")
        for ln in out_entries:
            if ln.strip().lower().startswith("##gff-version"):
                continue
            outf.write(ln.rstrip("\n") + "\n")

if __name__ == "__main__":
    main()
'''

# ---------------- fast_group_and_renumber_rename_feats.py (the second script you gave) ----------------
FAST_GROUP_PY = r'''#!/usr/bin/env python3
# fast_group_and_renumber_rename_feats.py
# (content identical to second script you provided)
from __future__ import annotations
import os, sys, re, argparse, tempfile, shutil
from collections import defaultdict, namedtuple
from pathlib import Path
import multiprocessing as mp

Gene = namedtuple("Gene", ["start","end","id","line"])

def parse_attrs(attrstr: str):
    d = {}
    if not attrstr or attrstr == ".":
        return d
    for part in attrstr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k,v = part.split("=",1)
            d[k] = v
        else:
            d[part] = ""
    return d

def build_attr(d: dict):
    if not d: return "."
    parts = []
    for k,v in d.items():
        if v == "": parts.append(k)
        else: parts.append(f"{k}={v}")
    return ";".join(parts)

def nat_seq_key(s):
    parts = re.split(r'(\d+)', s)
    key=[]
    for p in parts:
        if p.isdigit(): key.append((0,int(p)))
        else: key.append((1,p))
    return key

def find_gene_for_feature(gene_list, fstart, fend):
    if fstart is None or fend is None or not gene_list:
        return None
    lo = 0; hi = len(gene_list)-1; idx = None
    while lo <= hi:
        mid = (lo+hi)//2
        if gene_list[mid].start <= fstart:
            idx = mid; lo = mid + 1
        else:
            hi = mid - 1
    start_scan = max(0, (idx or 0) - 5)
    for i in range(start_scan, len(gene_list)):
        g = gene_list[i]
        if g.start is None:
            continue
        if g.start > fend:
            break
        if g.start <= fstart and g.end >= fend:
            return i
    return None

def sanitize_filename(s):
    return re.sub(r'[^A-Za-z0-9_.-]', '_', s)[:200]

def process_seq_bucket(args):
    seq, genes_sorted, tmpdir, out_path = args
    seq_tmp = os.path.join(tmpdir, f"out_{sanitize_filename(seq)}.gff")
    with open(seq_tmp, "w", encoding="utf-8") as outf:
        for g in genes_sorted:
            gid = g.id
            gene_file = os.path.join(tmpdir, f"{sanitize_filename(seq)}__{sanitize_filename(gid)}.gene")
            if os.path.exists(gene_file):
                with open(gene_file, "r", encoding="utf-8") as gf:
                    gene_line = gf.read().rstrip("\n")
            else:
                gene_line = f"{seq}\t.\tgene\t{g.start}\t{g.end}\t.\t.\t.\tID={gid};Name={gid}"
            outf.write(gene_line + "\n")

            mrna_file = os.path.join(tmpdir, f"{sanitize_filename(seq)}__{sanitize_filename(gid)}.mrna")
            mrnas = []
            if os.path.exists(mrna_file):
                with open(mrna_file, "r", encoding="utf-8") as mf:
                    for ln in mf:
                        line = ln.rstrip("\n")
                        cols = line.split("\t")
                        while len(cols) < 9: cols += ['.']
                        try:
                            s = int(cols[3])
                        except:
                            s = None
                        attrs = parse_attrs(cols[8])
                        oldid = attrs.get("ID") or attrs.get("id") or None
                        mrnas.append({"start": s, "oldid": oldid, "raw": line, "cols": cols, "attrs": attrs})
            mrnas.sort(key=lambda x: (x["start"] if x["start"] is not None else 10**12))

            base = gid
            if base.startswith("gene-"):
                base = base[len("gene-"):]
            enumerated = []
            old_to_new = {}
            for i, m in enumerate(mrnas, start=1):
                newid = f"{base}.t{int(i):02d}"
                if m["oldid"]:
                    old_to_new[m["oldid"]] = newid
                enumerated.append((newid, m))

            feat_file = os.path.join(tmpdir, f"{sanitize_filename(seq)}__{sanitize_filename(gid)}.feat")
            feats = []
            parent_to_feat_idxs = defaultdict(list)
            if os.path.exists(feat_file):
                with open(feat_file, "r", encoding="utf-8") as ff:
                    for j, ln in enumerate(ff):
                        ln = ln.rstrip("\n")
                        cols = ln.split("\t")
                        while len(cols) < 9: cols += ['.']
                        attrs = parse_attrs(cols[8])
                        pstr = attrs.get("Parent") or attrs.get("parent") or ""
                        parents = [p for p in pstr.split(",") if p] if pstr else []
                        feats.append({"raw": ln, "cols": cols, "attrs": attrs, "parents": parents})
                        for p in parents:
                            parent_to_feat_idxs[p].append(j)

            emitted_feats = set()

            if not enumerated:
                newid = f"{base}.t01"
                synth_attrs = {"ID": newid, "Parent": gid}
                synth_line = f"{seq}\t.\tmRNA\t{g.start}\t{g.end}\t.\t.\t.\t{build_attr(synth_attrs)}"
                outf.write(synth_line + "\n")
                counters = defaultdict(int)
                for fi, f in enumerate(feats):
                    cols_f = list(f["cols"])
                    attrs_f = dict(f["attrs"])
                    ftype = cols_f[2] if len(cols_f) > 2 else "feat"
                    counters[ftype] += 1
                    new_feat_id = f"{newid}.{ftype}{counters[ftype]}"
                    attrs_f['ID'] = new_feat_id
                    for pkey in ('Parent','parent'):
                        if pkey in attrs_f:
                            parts = attrs_f[pkey].split(',')
                            newparts = []
                            for p in parts:
                                if p in old_to_new:
                                    newparts.append(old_to_new[p])
                                elif p == gid:
                                    newparts.append(newid)
                                else:
                                    newparts.append(p)
                            attrs_f[pkey] = ",".join(newparts)
                            break
                    cols_f[8] = build_attr(attrs_f)
                    outf.write("\t".join(cols_f) + "\n")
                continue

            per_transcript_counters = defaultdict(lambda: defaultdict(int))

            for idx_enum, (newid, m) in enumerate(enumerated):
                cols = list(m["cols"])
                attrs = dict(m["attrs"])
                attrs['ID'] = newid
                attrs['Parent'] = gid
                cols[8] = build_attr(attrs)
                if cols[3] and cols[3] != '.': cols[3] = str(cols[3])
                if cols[4] and cols[4] != '.': cols[4] = str(cols[4])
                outf.write("\t".join(cols) + "\n")

                for fi, f in enumerate(feats):
                    if fi in emitted_feats:
                        continue
                    matched = False
                    for p in f["parents"]:
                        if m["oldid"] and p == m["oldid"]:
                            matched = True; break
                        if p in old_to_new and old_to_new[p] == newid:
                            matched = True; break
                        if p == gid and idx_enum == 0:
                            matched = True; break
                    if not matched:
                        continue
                    cols_f = list(f["cols"])
                    attrs_f = dict(f["attrs"])
                    ftype = cols_f[2] if len(cols_f) > 2 else "feat"
                    per_transcript_counters[newid][ftype] += 1
                    new_feat_id = f"{newid}.{ftype}{per_transcript_counters[newid][ftype]}"
                    attrs_f['ID'] = new_feat_id
                    for pkey in ('Parent','parent'):
                        if pkey in attrs_f:
                            parts = attrs_f[pkey].split(',')
                            newparts = []
                            for p in parts:
                                if p in old_to_new:
                                    newparts.append(old_to_new[p])
                                elif p == gid:
                                    newparts.append(enumerated[0][0])
                                else:
                                    newparts.append(p)
                            attrs_f[pkey] = ",".join(newparts)
                            break
                    cols_f[8] = build_attr(attrs_f)
                    outf.write("\t".join(cols_f) + "\n")
                    emitted_feats.add(fi)

            first_newid = enumerated[0][0] if enumerated else f"{base}.t01"
            for fi, f in enumerate(feats):
                if fi in emitted_feats:
                    continue
                cols_f = list(f["cols"])
                attrs_f = dict(f["attrs"])
                ftype = cols_f[2] if len(cols_f) > 2 else "feat"
                per_transcript_counters[first_newid][ftype] += 1
                new_feat_id = f"{first_newid}.{ftype}{per_transcript_counters[first_newid][ftype]}"
                attrs_f['ID'] = new_feat_id
                for pkey in ('Parent','parent'):
                    if pkey in attrs_f:
                        parts = attrs_f[pkey].split(',')
                        newparts = []
                        for p in parts:
                            if p in old_to_new:
                                newparts.append(old_to_new[p])
                            elif p == gid:
                                newparts.append(first_newid)
                            else:
                                newparts.append(p)
                        attrs_f[pkey] = ",".join(newparts)
                        break
                cols_f[8] = build_attr(attrs_f)
                outf.write("\t".join(cols_f) + "\n")

    return seq_tmp

def main():
    p = argparse.ArgumentParser(description="Fast group-by-gene renumbering with feature rename.")
    p.add_argument("in_gff", help="input gff3 (plain text, not gz)")
    p.add_argument("-o","--out", required=True, help="output gff3")
    p.add_argument("--tmpdir", default=None, help="temp directory")
    p.add_argument("--workers", type=int, default=1, help="parallel seq workers")
    args = p.parse_args()

    inp = args.in_gff
    out = args.out
    tmpdir = args.tmpdir or tempfile.mkdtemp(prefix="fastgff_")
    os.makedirs(tmpdir, exist_ok=True)
    workers = max(1, args.workers)

    genes_by_seq = defaultdict(list)
    header_lines = []
    print("PASS1: scanning genes ...", file=sys.stderr)
    with open(inp, "r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.strip():
                continue
            if ln.startswith("#"):
                if ln.strip().lower().startswith("##gff-version"):
                    header_lines.append("##gff-version 3")
                else:
                    header_lines.append(ln.rstrip("\n"))
                continue
            cols = ln.rstrip("\n").split("\t")
            while len(cols) < 9: cols += ['.']
            seqid = cols[0]; ftype = cols[2]
            try:
                s = int(cols[3]); e = int(cols[4])
            except:
                s = None; e = None
            attrs = parse_attrs(cols[8])
            if ftype.lower() == "gene":
                gid = attrs.get("ID") or attrs.get("id")
                if not gid:
                    gid = f"gene_{seqid}_{s}_{e}"
                    cols[8] = f"ID={gid};Name={gid}"
                g = Gene(start=s, end=e, id=gid, line="\t".join(cols))
                genes_by_seq[seqid].append(g)
                gfpath = os.path.join(tmpdir, f"{sanitize_filename(seqid)}__{sanitize_filename(gid)}.gene")
                with open(gfpath, "w", encoding="utf-8") as gf:
                    gf.write("\t".join(cols) + "\n")

    for seq in genes_by_seq:
        genes_by_seq[seq].sort(key=lambda x: (x.start if x.start is not None else -1))

    print("PASS2: dispatching features to gene buckets ...", file=sys.stderr)
    with open(inp, "r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            cols = ln.rstrip("\n").split("\t")
            while len(cols) < 9: cols += ['.']
            seqid = cols[0]; ftype = cols[2]
            try:
                s = int(cols[3]); e = int(cols[4])
            except:
                s = None; e = None
            attrs = parse_attrs(cols[8])
            if ftype.lower() == "gene":
                continue
            if seqid not in genes_by_seq or not genes_by_seq[seqid]:
                with open(os.path.join(tmpdir, "orphans.feat"), "a", encoding="utf-8") as of:
                    of.write(ln)
                continue
            idx = find_gene_for_feature(genes_by_seq[seqid], s, e)
            if idx is None:
                with open(os.path.join(tmpdir, "orphans.feat"), "a", encoding="utf-8") as of:
                    of.write(ln)
            else:
                gid = genes_by_seq[seqid][idx].id
                if ftype.lower() == "mrna":
                    with open(os.path.join(tmpdir, f"{sanitize_filename(seqid)}__{sanitize_filename(gid)}.mrna"), "a", encoding="utf-8") as mf:
                        mf.write(ln)
                else:
                    with open(os.path.join(tmpdir, f"{sanitize_filename(seqid)}__{sanitize_filename(gid)}.feat"), "a", encoding="utf-8") as ff:
                        ff.write(ln)

    print("FINALIZING per-seq buckets ...", file=sys.stderr)
    seqs = sorted(genes_by_seq.keys(), key=lambda x: nat_seq_key(x))
    tasks = [(seq, genes_by_seq[seq], tmpdir, out) for seq in seqs]

    seq_out_files = []
    if workers > 1:
        with mp.Pool(processes=workers) as pool:
            for res in pool.imap_unordered(process_seq_bucket, tasks):
                seq_out_files.append(res)
    else:
        for t in tasks:
            seq_out_files.append(process_seq_bucket(t))

    print("WRITING final output ...", file=sys.stderr)
    with open(out, "w", encoding="utf-8") as outf:
        seen_ver = False
        for h in header_lines:
            if h.strip().lower().startswith("##gff-version"):
                if not seen_ver:
                    outf.write("##gff-version 3\n"); seen_ver = True
            else:
                outf.write(h.rstrip("\n") + "\n")
        for seq in seqs:
            seq_tmp = os.path.join(tmpdir, f"out_{sanitize_filename(seq)}.gff")
            if os.path.exists(seq_tmp):
                with open(seq_tmp, "r", encoding="utf-8") as st:
                    shutil.copyfileobj(st, outf)
        orp = os.path.join(tmpdir, "orphans.feat")
        if os.path.exists(orp):
            outf.write("## other_orphan_features\n")
            with open(orp, "r", encoding="utf-8") as of:
                shutil.copyfileobj(of, outf)

    print("Done. output ->", out, file=sys.stderr)
    # shutil.rmtree(tmpdir)  # uncomment to auto-clean temp files

if __name__ == "__main__":
    main()
'''

# ---------------- prepare and submit per-genome (main) ----------------
def prepare_and_submit(genome: str, protein: str, input_dir: str,
                       partition: str, ntasks_per_node: int, threads: int, job_name_prefix: str, group_workers: int):
    genome = os.path.abspath(genome)
    protein = os.path.abspath(protein)
    prefix = Path(genome).stem
    job_dir = os.path.join(input_dir, prefix)
    safe_mkdir(job_dir)

    index_file = os.path.join(job_dir, "genome.mpi")
    raw_gff = os.path.join(job_dir, "aln.gff")
    fixed_step2 = os.path.join(job_dir, "aln.fixed.step2.gff")
    strict_gff = os.path.join(job_dir, "aln.strict.gff")
    final_gff = os.path.join(input_dir, f"{prefix}.gff3")
    sbatch_path = os.path.join(job_dir, "run_miniprot.sbatch")
    fix_script_path = os.path.join(job_dir, "fix_and_ensure.py")
    strict_script_path = os.path.join(job_dir, "strict_filter_merge.py")
    group_script_path = os.path.join(job_dir, "fast_group_and_renumber_rename_feats.py")
    job_name = job_name_prefix if job_name_prefix else f"pep_{prefix}"

    # write helper scripts
    with open(fix_script_path, "w", encoding="utf-8") as fh:
        fh.write(FIX_AND_ENSURE_PY)
    with open(strict_script_path, "w", encoding="utf-8") as fh:
        fh.write(STRICT_FILTER_MERGE_PY)
    with open(group_script_path, "w", encoding="utf-8") as fh:
        fh.write(FAST_GROUP_PY)

    os.chmod(fix_script_path, 0o750)
    os.chmod(strict_script_path, 0o750)
    os.chmod(group_script_path, 0o750)
    print(f"[WROTE] helper scripts to {job_dir}")

    content = SBATCH_TEMPLATE.format(
        job_name=job_name,
        partition=partition,
        ntasks_per_node=ntasks_per_node,
        job_dir=job_dir,
        genome=genome,
        protein=protein,
        threads=threads,
        index_file=index_file,
        raw_gff=raw_gff,
        fixed_step2=fixed_step2,
        strict_gff=strict_gff,
        final_gff=final_gff,
        fix_script=fix_script_path,
        strict_script=strict_script_path,
        group_script=group_script_path,
        group_workers=group_workers
    )

    with open(sbatch_path, "w", encoding="utf-8") as fh:
        fh.write(content)
    os.chmod(sbatch_path, 0o750)
    print(f"[WROTE] {sbatch_path}")

    # Submit immediately
    try:
        res = subprocess.run(["sbatch", sbatch_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out = res.stdout.strip()
        print(f"[SUBMITTED] {prefix}: {out}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] sbatch failed for {prefix}: {e.stderr}", file=sys.stderr)

# ---------------- main (CLI) ----------------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Batch-submit miniprot + strict_filter_merge + fast_group jobs.")
    parser.add_argument("--input-dir", required=True, help="Directory containing genome fasta files (*.fna/*.fa/*.fasta)")
    parser.add_argument("--protein", required=True, help="Protein fasta (miniprotzhushi.fa)")
    parser.add_argument("--partition", default="hebhcnormal01", help="SBATCH partition (default hebhcnormal01)")
    parser.add_argument("--ntasks-per-node", type=int, default=60, help="SBATCH --ntasks-per-node (default 60)")
    parser.add_argument("--threads", type=int, default=16, help="Threads for miniprot (default 16)")
    parser.add_argument("--job-name", default="pep", help="SBATCH job-name (default pep)")
    parser.add_argument("--group-workers", type=int, default=1, help="workers passed to grouping script (per-job parallelism)")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    protein = os.path.abspath(args.protein)
    partition = args.partition
    ntasks = args.ntasks_per_node
    threads = args.threads
    jobname = args.job_name
    group_workers = max(1, args.group_workers)

    if not os.path.isdir(input_dir):
        print("Input directory not found:", input_dir, file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(protein):
        print("Protein fasta not found:", protein, file=sys.stderr)
        sys.exit(1)

    genome_files = find_genome_files(input_dir)
    if not genome_files:
        print("No genome fasta files found in:", input_dir)
        sys.exit(1)

    print(f"Found {len(genome_files)} genomes in {input_dir}. Preparing jobs...")
    for g in genome_files:
        prepare_and_submit(genome=g, protein=protein, input_dir=input_dir,
                           partition=partition, ntasks_per_node=ntasks,
                           threads=threads, job_name_prefix=jobname, group_workers=group_workers)

    print("All jobs prepared and submitted.")

if __name__ == "__main__":
    main()

