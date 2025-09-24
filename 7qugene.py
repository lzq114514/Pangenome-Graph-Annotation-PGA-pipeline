#!/usr/bin/env python3
"""
strict_filter_merge.py

严格实现三条规则（identical / partial-overlap delete / containment -> reparent as transcripts）。
Two-pass: parse -> decide -> output (safe, avoids index-shift bugs).
支持 plain gff3 和 .gff3.gz。
"""

import sys
import gzip
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
        elif " " in part:
            k,v = part.split(" ",1)
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

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 strict_filter_merge.py input.gff3[.gz] > out.gff3", file=sys.stderr)
        sys.exit(1)
    path = sys.argv[1]

    # 1) read all lines into list of dicts (preserve original order)
    entries = []
    try:
        fh = open_maybe_gz(path)
    except Exception as e:
        print("Error opening file: {}".format(e), file=sys.stderr)
        sys.exit(2)

    with fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#") or not line.strip():
                entries.append({"is_feature": False, "raw": line})
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                entries.append({"is_feature": False, "raw": line})
                continue
            seqid = cols[0]
            ftype = cols[2]
            try:
                start = int(cols[3])
                end = int(cols[4])
            except:
                start = None
                end = None
            attrstr = cols[8] if len(cols) > 8 else "."
            attrs = parse_attrs(attrstr)
            # normalize ID/Parent keys retrieval
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

    # 2) Build id -> list(indices) and parent -> children indices
    id_to_idxs = defaultdict(list)
    parent_to_children = defaultdict(list)
    for i,e in enumerate(entries):
        if not e.get("is_feature"):
            continue
        if e.get("id"):
            id_to_idxs[e["id"]].append(i)
        for p in e.get("parents", []):
            parent_to_children[p].append(i)

    # 3) Collect gene list per seqid in file order
    genes_by_seq = defaultdict(list)  # seq -> list of tuples (idx, id, start, end)
    for i,e in enumerate(entries):
        if not e.get("is_feature"): continue
        if e["type"].lower() == "gene":
            genes_by_seq[e["seqid"]].append((i, e.get("id"), e.get("start"), e.get("end")))

    # 4) Decide: which genes to keep, which to remove, which are contained->map child->parent
    kept_genes_by_seq = defaultdict(list)  # seq -> list of (idx,id,s,e) for kept genes in file order
    removed_gene_ids = set()               # genes to DELETE with all descendants (identical or overlap)
    contained_map = {}                     # child_gene_id -> parent_gene_id (strict containment)
    # We'll also ensure that identical-case and overlap-case result in deletion of later gene and its descendants

    for seq, glist in genes_by_seq.items():
        for idx,gid,s,e in glist:
            # if missing info, conservatively keep
            if gid is None or s is None or e is None:
                kept_genes_by_seq[seq].append((idx,gid,s,e))
                continue
            # compare against all previously kept genes on this seq
            identical = False
            contained_by = None
            overlaps_any = False
            for ks, kid, kstart, kend in [(x[0], x[1], x[2], x[3]) for x in kept_genes_by_seq[seq]]:
                # note: we stored kept as (idx,id,s,e) but iterate in same order; get values properly:
                pass
            # better: iterate stored list
            for (kidx, kid, ks, kt) in kept_genes_by_seq[seq]:
                if ks is None or kt is None:
                    continue
                if s == ks and e == kt:
                    identical = True
                    break
                if is_strictly_contained(s,e,ks,kt):
                    # child completely inside kept gene -> contained
                    contained_by = kid
                    break
                if intervals_overlap(s,e,ks,kt):
                    overlaps_any = True
            if identical:
                removed_gene_ids.add(gid)
                continue
            if contained_by:
                # map child -> parent (do NOT add to kept)
                contained_map[gid] = contained_by
                continue
            if overlaps_any:
                removed_gene_ids.add(gid)
                continue
            # otherwise keep
            kept_genes_by_seq[seq].append((idx,gid,s,e))

    # 5) For contained genes, check whether they have mRNA children.
    #    If not, but they have direct children (CDS/exon etc), we will plan to synthesize one synthetic mRNA per child-gene.
    need_synth_mrna_for = {}  # child_gene_id -> synthetic_mrna_id (to be generated lazily)
    synth_counter = 0
    def make_synth_id(parent_gene_id):
        nonlocal synth_counter
        synth_counter += 1
        return f"{parent_gene_id}.auto_mrna{synth_counter}"

    for child_gid, parent_gid in contained_map.items():
        # find children indices whose parent include child_gid
        child_idxs = parent_to_children.get(child_gid, [])
        has_mrna = False
        has_other = False
        for ci in child_idxs:
            ce = entries[ci]
            if ce["type"].lower() == "mrna":
                has_mrna = True
            else:
                # skip if it's the gene line itself
                if ce["type"].lower() == "gene":
                    continue
                has_other = True
        if not has_mrna and has_other:
            # need to plan a synthetic mRNA to hold these other features
            need_synth_mrna_for[child_gid] = None  # id allocated later on-demand

    # 6) Build set of all entry indices to remove because they descend from removed_gene_ids
    #    BFS on parent->children using gene ids as starting points
    to_remove_entry_idxs = set()
    # seed with gene entry indices whose id in removed_gene_ids
    q = deque()
    for gid in removed_gene_ids:
        for gi in id_to_idxs.get(gid, []):
            # only gene entries
            ge = entries[gi]
            if ge["type"].lower() == "gene":
                q.append(gi)
                to_remove_entry_idxs.add(gi)
    # BFS collect descendants
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

    # Additionally: we must delete the child gene lines themselves when they are in contained_map (they are not "removed" as overlap/identical,
    # but per rule 3 child gene lines should NOT be written; so mark their gene-entry idx for removal)
    for child_gid in contained_map.keys():
        for gi in id_to_idxs.get(child_gid, []):
            # if gene line, mark removal (but will reparent its children instead)
            ge = entries[gi]
            if ge["type"].lower() == "gene":
                to_remove_entry_idxs.add(gi)

    # 7) Now do final output pass: iterate original entries in order and write according to rules:
    #    - comments/empty -> print
    #    - if entry idx in to_remove_entry_idxs -> skip
    #    - if entry's Parent includes a contained child gene id -> reparent:
    #         a) if parent gene was 'contained' (we have contained_map[child]->parent), then:
    #             - if entry is mRNA: change its Parent from child_gid to parent_gid
    #             - else (e.g., CDS/exon directly under child): if synthetic mRNA exists for child_gid, reparent to synthetic mRNA; else create synthetic mRNA now, emit it, record it, then reparent
    #    - otherwise emit entry (but ensure attr string updated if Parent changed)
    #
    # We'll need to know whether we've already emitted a synthetic mRNA for particular child_gid (so not to emit twice).
    emitted_synth = {}  # child_gid -> synthetic_mrna_id (once created)
    parent_map = contained_map  # alias

    out_lines = []
    for idx,e in enumerate(entries):
        if not e.get("is_feature"):
            out_lines.append(e["raw"])
            continue
        # skip entries marked for removal
        if idx in to_remove_entry_idxs:
            continue
        # process features whose parents include child_gene ids (that were contained)
        cur_parents = list(e["parents"]) if e.get("parents") else []
        changed = False
        new_parents = []
        for p in cur_parents:
            if p in parent_map:
                # this feature originally parented to a contained child gene
                parent_gene = parent_map[p]  # the gene that will keep the transcripts
                if e["type"].lower() == "mrna":
                    # reparent mRNA directly to parent_gene
                    new_parents.append(parent_gene)
                    changed = True
                else:
                    # non-mRNA child: need a synthetic mRNA under parent_gene
                    # create (or reuse) synthetic mRNA id for this child gene
                    if emitted_synth.get(p) is None:
                        syn_id = make_synth_id(parent_gene)
                        emitted_synth[p] = syn_id
                        # create a synthetic mRNA line and emit it immediately BEFORE this entry
                        # derive coords from this feature and possibly child siblings (we'll use this feature's coords as mRNA coords)
                        s = e.get("start") or entries[id_to_idxs[p][0]]["start"] if id_to_idxs.get(p) else 0
                        t = e.get("end")   or entries[id_to_idxs[p][0]]["end"]   if id_to_idxs.get(p) else 0
                        syn_attrs = {"ID": syn_id, "Parent": parent_gene}
                        syn_attrstr = build_attrstr(syn_attrs)
                        syn_cols = [e["seqid"], ".", "mRNA", str(s), str(t), ".", ".", ".", syn_attrstr]
                        out_lines.append("\t".join(syn_cols))
                    # reparent this feature to the synthetic mRNA id
                    new_parents.append(emitted_synth[p])
                    changed = True
            else:
                new_parents.append(p)
        if changed:
            # update attrs: replace Parent (prefer original key capitalization if present)
            # prefer to update existing 'Parent' or 'parent' key; else add 'Parent'
            attrs = dict(e["attrs"])  # copy
            if "Parent" in attrs:
                attrs["Parent"] = ",".join(new_parents)
            elif "parent" in attrs:
                attrs["parent"] = ",".join(new_parents)
            else:
                attrs["Parent"] = ",".join(new_parents)
            # rebuild attr string and raw line
            cols = list(e["cols"])
            if len(cols) > 8:
                cols[8] = build_attrstr(attrs)
            else:
                while len(cols) < 8:
                    cols.append(".")
                cols.append(build_attrstr(attrs))
            out_lines.append("\t".join(cols))
        else:
            # normal emit (raw)
            out_lines.append(e["raw"])

    # write all out
    sys.stdout.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
