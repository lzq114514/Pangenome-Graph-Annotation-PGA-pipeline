#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_miniprot_submit_all.py

Batch-submit miniprot jobs for all genome fasta files in a directory.
Each genome -> a subfolder under out_dir -> a sbatch script is written AND SUBMITTED.

Defaults (as requested):
#SBATCH --job-name=pep
#SBATCH --partition=hebhcnormal01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --error=%j.err
#SBATCH --output=%j.out

Usage example:
  python3 batch_miniprot_submit_all.py \
    --input-dir /path/to/genomes \
    --protein /path/to/miniprotzhushi.fa \
    --out-dir /path/to/out \
    --cpus-per-task 16

Notes:
 - This script WILL submit every generated sbatch script (no --no-submit option).
 - It writes per-sample job dirs under out_dir/<genome-stem>/ with run_miniprot.sbatch.
 - Each job builds miniprot index, runs miniprot to produce aln.gff, then runs an inline
   Python "fix_gff" routine to create pga_anno.gff (cleaned GFF3).
 - Index, raw_gff and fixed_gff are overwritten if they exist.
"""

from __future__ import annotations
import os
import sys
import glob
import subprocess
from pathlib import Path
from typing import List

# ---------------- utility functions ----------------
def find_genome_files(input_dir: str) -> List[str]:
    exts = ("*.fna", "*.fa", "*.fasta", "*.FNA", "*.FA", "*.FASTA")
    files = []
    for e in exts:
        files.extend(sorted(glob.glob(os.path.join(input_dir, e))))
    return files

def safe_mkdir(path: str):
    os.makedirs(path, exist_ok=True)

# ---------------- SBATCH template (user-specified defaults) ----------------
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
echo "Fixed GFF: {fixed_gff}"

# Build index
miniprot -t{threads} -d "{index_file}" "{genome}"

# Align and write raw GFF
miniprot -t{threads} --gff "{index_file}" "{protein}" > "{raw_gff}"

# Inline GFF fixer (creates pga_anno.gff)
python3 - <<'PY'
import re, sys
in_gff = r"{raw_gff}"
out_gff = r"{fixed_gff}"

def fix_gff(in_path, out_path):
    gene_count = 1
    exon_count = 1
    current_parent = None
    with open(in_path, 'r', encoding='utf-8', errors='replace') as fin, open(out_path, 'w', encoding='utf-8') as fout:
        print("##gff-version 3", file=fout)
        for line in fin:
            if not line.strip() or line.startswith("#") or line.startswith("##PAF"):
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attr = cols[:9]

            if feature == "mRNA":
                exon_count = 1
                m = re.search(r'ID=([^;\\s]+)', attr)
                if m:
                    mrna_id = m.group(1)
                else:
                    mrna_id = f"MP{{gene_count:06d}}"
                    attr = f"ID={{mrna_id}};" + attr
                gene_id = f"gene-{{mrna_id}}"
                gene_attr = f"ID={{gene_id}};Name={{gene_id}}"
                print("\\t".join([seqid, source, "gene", start, end, score, strand, ".", gene_attr]), file=fout)
                if 'Parent=' in attr:
                    attr = re.sub(r'Parent=[^;\\s]+', f"Parent={{gene_id}}", attr)
                else:
                    if 'ID=' not in attr:
                        attr = f"ID={{mrna_id}};Parent={{gene_id}};" + attr
                    else:
                        attr = re.sub(r'ID=[^;\\s]+', f"ID={{mrna_id}};Parent={{gene_id}}", attr)
                print("\\t".join([seqid, source, "mRNA", start, end, score, strand, ".", attr]), file=fout)
                current_parent = mrna_id
                gene_count += 1
                continue

            if feature == "exon" and current_parent:
                if 'ID=' not in attr:
                    exon_id = f"exon-{{current_parent}}-{{exon_count}}"
                    exon_count += 1
                    attr = f"ID={{exon_id}};Parent={{current_parent}}"
                else:
                    if 'Parent=' in attr:
                        attr = re.sub(r'Parent=[^;\\s]+', f"Parent={{current_parent}}", attr)
                    else:
                        attr = attr + f";Parent={{current_parent}}"
                print("\\t".join([seqid, source, "exon", start, end, score, strand, phase, attr]), file=fout)
                continue

            if feature == "CDS" and current_parent:
                exon_id = f"exon-{{current_parent}}-{{exon_count}}"
                exon_count += 1
                exon_attr = f"ID={{exon_id}};Parent={{current_parent}}"
                print("\\t".join([seqid, source, "exon", start, end, ".", strand, ".", exon_attr]), file=fout)

                if 'Parent=' in attr:
                    attr_cds = re.sub(r'Parent=[^;\\s]+', f"Parent={{current_parent}}", attr)
                else:
                    attr_cds = f"Parent={{current_parent}};" + attr
                print("\\t".join([seqid, source, "CDS", start, end, score, strand, phase, attr_cds]), file=fout)
                continue

            if feature == "stop_codon" and current_parent:
                if 'Parent=' in attr:
                    attr_stop = re.sub(r'Parent=[^;\\s]+', f"Parent={{current_parent}}", attr)
                else:
                    attr_stop = f"Parent={{current_parent}};" + attr
                print("\\t".join([seqid, source, "stop_codon", start, end, score, strand, phase, attr_stop]), file=fout)
                continue

            print("\\t".join([seqid, source, feature, start, end, score, strand, phase, attr]), file=fout)

try:
    fix_gff(in_gff, out_gff)
    print("GFF fixed:", out_gff)
except Exception as e:
    print("GFF fix failed:", e, file=sys.stderr)
    sys.exit(2)
PY

echo "Job done: $(date)"
"""

# ---------------- prepare and submit per-genome ----------------
def prepare_and_submit(genome: str, protein: str, out_dir: str,
                       partition: str, ntasks_per_node: int, cpus_per_task: int,
                       threads: int, job_name_prefix: str, extra_sbatch: str):
    genome = os.path.abspath(genome)
    protein = os.path.abspath(protein)
    prefix = Path(genome).stem
    job_dir = os.path.join(out_dir, prefix)
    safe_mkdir(job_dir)

    index_file = os.path.join(job_dir, "genome.mpi")
    raw_gff = os.path.join(job_dir, "aln.gff")
    fixed_gff = os.path.join(job_dir, "pga_anno.gff")
    sbatch_path = os.path.join(job_dir, "run_miniprot.sbatch")
    job_name = job_name_prefix if job_name_prefix else f"pep_{prefix}"

    extra = extra_sbatch or ""

    content = SBATCH_TEMPLATE.format(
        job_name=job_name,
        partition=partition,
        ntasks_per_node=ntasks_per_node,
        extra_sbatch=extra,
        job_dir=job_dir,
        genome=genome,
        protein=protein,
        threads=threads,
        index_file=index_file,
        raw_gff=raw_gff,
        fixed_gff=fixed_gff
    )

    with open(sbatch_path, "w", encoding="utf-8") as fh:
        fh.write(content)
    os.chmod(sbatch_path, 0o750)
    print(f"[WROTE] {sbatch_path}")

    # Submit immediately (per user's request defaults to submit all)
    try:
        res = subprocess.run(["sbatch", sbatch_path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out = res.stdout.strip()
        print(f"[SUBMITTED] {prefix}: {out}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] sbatch failed for {prefix}: {e.stderr}", file=sys.stderr)

# ---------------- main (CLI) ----------------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Batch-submit miniprot jobs (auto-submit).")
    parser.add_argument("--input-dir", required=True, help="Directory containing genome fasta files (*.fna/*.fa/*.fasta)")
    parser.add_argument("--protein", required=True, help="Protein fasta (miniprotzhushi.fa)")
    parser.add_argument("--out-dir", default="./miniprot_batch_out", help="Base output directory for jobs")
    parser.add_argument("--partition", default="hebhcnormal01", help="SBATCH partition (default hebhcnormal01)")
    parser.add_argument("--ntasks-per-node", type=int, default=60, help="SBATCH --ntasks-per-node (default 60)")
    parser.add_argument("--cpus-per-task", type=int, default=16, help="SBATCH --cpus-per-task (default 16)")
    parser.add_argument("--threads", type=int, default=None, help="Threads for miniprot (default = cpus-per-task)")
    parser.add_argument("--job-name", default="pep", help="SBATCH job-name (default pep)")
    parser.add_argument("--extra-sbatch", default="", help="Extra SBATCH lines (raw text, optional)")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    protein = os.path.abspath(args.protein)
    out_dir = os.path.abspath(args.out_dir)
    partition = args.partition
    ntasks = args.ntasks_per_node
    cpus = args.cpus_per_task
    threads = args.threads if args.threads is not None else cpus
    jobname = args.job_name
    extra = args.extra_sbatch

    if not os.path.isdir(input_dir):
        print("Input directory not found:", input_dir, file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(protein):
        print("Protein fasta not found:", protein, file=sys.stderr)
        sys.exit(1)
    safe_mkdir(out_dir)

    genome_files = find_genome_files(input_dir)
    if not genome_files:
        print("No genome fasta files found in:", input_dir)
        sys.exit(1)

    print(f"Found {len(genome_files)} genomes. Preparing jobs in: {out_dir}")
    for g in genome_files:
        prepare_and_submit(genome=g, protein=protein, out_dir=out_dir,
                           partition=partition, ntasks_per_node=ntasks,
                           cpus_per_task=cpus, threads=threads,
                           job_name_prefix=jobname, extra_sbatch=extra)

    print("All jobs prepared and submitted.")

if __name__ == "__main__":
    main()
