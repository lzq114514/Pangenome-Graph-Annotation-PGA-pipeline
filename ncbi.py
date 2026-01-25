#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
reorder_gff3_exon_cds.py

Same logic as before, but processes a directory (or a single file) and
overwrites each original .gff3 file in-place (atomic replace).
"""

import sys
from pathlib import Path
from collections import defaultdict
import tempfile
import shutil
import os

STOP_FEATURES = {"stop", "stop_codon"}


# ---------- basic utils ----------

def parse_attrs(s):
    d = {}
    for p in s.split(";"):
        if "=" in p:
            k, v = p.split("=", 1)
            d[k] = v
    return d


def format_attrs(d):
    out = []
    if "ID" in d:
        out.append(f"ID={d['ID']}")
    if "Parent" in d:
        out.append(f"Parent={d['Parent']}")
    for k, v in d.items():
        if k not in ("ID", "Parent"):
            out.append(f"{k}={v}")
    return ";".join(out)


def overlap_len(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def overlap_ratio(a, b):
    ov = overlap_len(a, b)
    ml = min(a[1] - a[0] + 1, b[1] - b[0] + 1)
    return 0.0 if ml == 0 else ov / ml


def parse_phase(cols):
    ph = cols[7]
    if ph in (".", ""):
        return None
    try:
        p = int(ph)
        return p if p in (0, 1, 2) else None
    except Exception:
        return None


# ---------- boundary ----------

def find_other_boundary(lines):
    for i, ln in enumerate(lines):
        if ln.startswith("#") and ln.lstrip("#").strip().lower().startswith("other"):
            return i
    for i, ln in enumerate(lines):
        if ln.startswith("#"):
            continue
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2].lower() == "other":
            return i
    return None


# ---------- mRNA dedup (SAFE) ----------

def precompute_mrna_drop(lines):
    """
    Only mark mRNA lines to drop.
    gene lines are NEVER touched here.
    """
    mrna_info = []
    for i, ln in enumerate(lines):
        if ln.startswith("#"):
            continue
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 9 and c[2].lower() == "mrna":
            mrna_info.append((i, c))

    # index → next non-comment feature
    idx2type = {}
    for i, ln in enumerate(lines):
        if ln.startswith("#"):
            idx2type[i] = None
            continue
        c = ln.rstrip("\n").split("\t")
        idx2type[i] = c[2].lower() if len(c) >= 3 else None

    seen = {}
    drop = set()

    for i, c in mrna_info:
        nxt = i + 1
        # immune if immediate next line in pre is exon
        if nxt < len(lines) and idx2type.get(nxt) == "exon":
            continue  # immune

        key = (c[0], c[3], c[4], c[6])
        if key in seen:
            drop.add(i)
        else:
            seen[key] = i

    return drop


# ---------- main processing ----------

def process_file(infile_path, outfile_path):
    # Read as text preserving line endings
    text = Path(infile_path).read_text()
    lines = text.splitlines(keepends=True)

    boundary = find_other_boundary(lines)
    pre = lines if boundary is None else lines[:boundary]
    post = [] if boundary is None else lines[boundary:]

    mrna_drop = precompute_mrna_drop(pre)

    out = []
    current_mrna = None
    exon_cds = defaultdict(lambda: {"exon": None, "cds": None})

    def flush():
        nonlocal exon_cds
        if not current_mrna:
            exon_cds.clear()
            return

        items = []
        for (s, e) in sorted(exon_cds):
            pair = exon_cds[(s, e)]
            items.append({
                "coord": (s, e),
                "exon": pair["exon"],
                "cds": pair["cds"],
                "phase": parse_phase(pair["cds"]) if pair["cds"] else None
            })

        # phase-aware CDS overlap filter
        kept = []
        for iv in items:
            bad = False
            if iv["cds"]:
                for pv in kept:
                    if pv["cds"] and overlap_len(iv["coord"], pv["coord"]) > 0:
                        if iv["phase"] is not None and pv["phase"] is not None:
                            if iv["phase"] != pv["phase"]:
                                bad = True
                                break
            if not bad:
                kept.append(iv)

        # overlap collapse (keep first on any overlap)
        collapsed = []
        for iv in kept:
            if not collapsed:
                collapsed.append(iv)
            elif overlap_ratio(collapsed[-1]["coord"], iv["coord"]) > 0:
                continue
            else:
                collapsed.append(iv)

        ei = ci = 0
        for iv in collapsed:
            if iv["exon"]:
                ei += 1
                c = iv["exon"]
                a = parse_attrs(c[8])
                a["ID"] = f"{current_mrna}.exon{ei}"
                a["Parent"] = current_mrna
                c[8] = format_attrs(a)
                out.append("\t".join(c) + "\n")
            if iv["cds"]:
                ci += 1
                c = iv["cds"]
                a = parse_attrs(c[8])
                a["ID"] = f"{current_mrna}.CDS{ci}"
                a["Parent"] = current_mrna
                c[8] = format_attrs(a)
                out.append("\t".join(c) + "\n")

        exon_cds.clear()

    i = 0
    while i < len(pre):
        ln = pre[i]
        if ln.startswith("#"):
            out.append(ln)
            i += 1
            continue

        c = ln.rstrip("\n").split("\t")
        if len(c) < 9:
            out.append(ln)
            i += 1
            continue

        f = c[2].lower()

        if f in STOP_FEATURES:
            i += 1
            continue

        if f == "gene":
            flush()
            current_mrna = None
            out.append(ln)
            i += 1
            continue

        if f == "mrna":
            flush()
            if i in mrna_drop:
                # skip the mRNA line (drop transcript entirely)
                current_mrna = None
                i += 1
                # skip its child lines until next gene/mRNA/other boundary
                while i < len(pre):
                    nxt = pre[i]
                    if nxt.startswith("#"):
                        i += 1
                        continue
                    nxt_cols = nxt.rstrip("\n").split("\t")
                    if len(nxt_cols) >= 3:
                        nxt_ft = nxt_cols[2].lower()
                        # stop skipping when we reach next mRNA or next gene or an 'other' data line
                        if nxt_ft in ("mrna", "gene", "other"):
                            break
                    i += 1
                continue
            current_mrna = parse_attrs(c[8]).get("ID")
            out.append(ln)
            i += 1
            continue

        if f in ("exon", "cds"):
            if not current_mrna:
                # orphan child (no active mRNA) — skip to avoid emitting orphans
                i += 1
                continue
            try:
                key = (int(c[3]), int(c[4]))
            except Exception:
                out.append(ln)
                i += 1
                continue
            if exon_cds[key][f] is None:
                exon_cds[key][f] = c
            i += 1
            continue

        out.append(ln)
        i += 1

    flush()
    out.extend(post)

    # write into temporary file then atomically replace original
    temp_dir = Path(outfile_path).parent
    fd, tmp_path = tempfile.mkstemp(dir=str(temp_dir), prefix=Path(outfile_path).name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as tmpf:
            tmpf.write("".join(out))
        # replace original atomically
        Path(tmp_path).replace(outfile_path)
    finally:
        # ensure no leftover temp file (replace above removes tmp), but be safe:
        if Path(tmp_path).exists():
            try:
                Path(tmp_path).unlink()
            except Exception:
                pass


# ---------- entry ----------

def process_path(path_str):
    p = Path(path_str)
    if p.is_dir():
        gffs = sorted(p.glob("*.gff3"))
        if not gffs:
            print(f"No .gff3 files found in directory: {p}", file=sys.stderr)
            return
        for f in gffs:
            print(f"Processing (in-place) {f} ...")
            process_file(str(f), str(f))
        print("Done.")
    else:
        f = Path(p)
        if not f.exists() or f.suffix.lower() != ".gff3":
            print("Input must be a directory or a .gff3 file", file=sys.stderr)
            sys.exit(1)
        print(f"Processing (in-place) {f} ...")
        process_file(str(f), str(f))
        print("Done.")


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 reorder_gff3_exon_cds.py <gff3 file or dir>")
    process_path(sys.argv[1])


if __name__ == "__main__":
    main()
