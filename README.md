# graph pangenome annotation pipeline (PAP)

"A pggb-based homologous backfill annotation method"

## Table of Contents
1. [Setting up the environment](#setting-up-the-environment)  
2. [Input requirements](#input-requirements)  
3. [Running the pipeline](#running-the-pipeline)  
   - [Step 1: Prepare and run pggb](#step-1-prepare-and-run-pggb)  
   - [Step 2: Run PAP](#step-2-run-PAP)  
---

## Setting up the environment

Create a conda environment from the provided `pap.yaml`:


`conda env create -f pap.yaml`
Description:
This sets up all dependencies required to run the pipeline.

## Input requirements
A folder containing all genome FASTA files (.fna) and corresponding GFF files.

Chromosome IDs in FASTA files must follow the format "chromosome 1", "chromosome X" for sex chromosomes, etc. (Currently only human genome fully supported; the example uses Arabidopsis in data/ninanjie/fna).

<img width="270" height="70" alt="image" src="https://github.com/user-attachments/assets/72292c0b-e1d5-489a-9525-e8c5439ec979" />

## Running the pipeline
## Step 1: Prepare and run pggb
Run:


`sbatch /public/home/acfurbn1nz/huitian/github/PAP/step1-2.sh /public/home/acfurbn1nz/your_input_dir -p hebhcnormal01 -r ninanjie1 -s 1 -e 5 -t 20 --pggb-s 5000 --pggb-l 25000 --pggb-p 90 --pggb-c 1 --pggb-K 19 --pggb-F 0.001 --pggb-g 30 --pggb-k 23 --pggb-f 0 --pggb-B 10M --pggb-n 5 --pggb-j 0 --pggb-e 0 --pggb-G 700,900,1100 --pggb-P 1,19,39,3,81,1 --pggb-O 0.001 --pggb-d 100 --pggb-Q Consensus_ --pggb-Y "#"`

Description:
Parameter Description

Required
input_dir
Directory containing all genome assemblies and corresponding annotation files.
Each subdirectory should include files named as:

genome_name/genome_name.fna
genome_name/genome_name.gff

Example:
/public/home/acfurbn1nz/huitian/github/shuidao

⚙️ General Options
-p, --partition
SLURM partition (queue) name used for job submission.
Example: hebhcnormal01
(default: hebhcnormal01)

-r, --ref
Reference genome prefix used in PGGB graph construction.
Typically formatted as:

<input_dir_name> + "1"

Example:
shuidao1

-s, --start
Starting chromosome index.
Usually begins from 1.
(default: 1)
-e, --end
Ending chromosome index.
Should match the total chromosome number of the species.
(default: 5)
-t, --threads
Number of threads allocated per PGGB job.
(default: 20)
🧬 PGGB Parameters

These parameters control the behavior of PGGB during graph construction and alignment.

📌 Core alignment & graph construction
--pggb-s
Segment length used for sequence partitioning.
Larger values increase speed but may reduce sensitivity.
--pggb-l
Minimum alignment length.
Filters out short alignments.
--pggb-p
Minimum sequence identity (%) required for alignment.
--pggb-c
Minimum coverage required for including an alignment.
🔬 Mapping & seeding
--pggb-K
k-mer size used for initial mapping (seed length).
--pggb-F
Frequency threshold for filtering repetitive k-mers.
--pggb-g
Maximum gap length allowed in alignments.
--pggb-k
k-mer size used in graph construction (refinement stage).
--pggb-f
Additional filtering parameter for alignment refinement.
⚙️ Graph building behavior
--pggb-B
Memory limit (e.g., 10M) for graph partitioning blocks.
--pggb-n
Number of haplotypes/genomes expected in the graph.
--pggb-j
Number of parallel jobs used internally by PGGB.
--pggb-e
Edge filtering threshold for graph simplification.
🧩 Multi-scale alignment tuning
--pggb-G
Multi-scale alignment segment sizes (comma-separated).
Example: 700,900,1100
--pggb-P
Multi-stage alignment parameters controlling sensitivity/speed trade-offs.
📊 Output & refinement
--pggb-O
Alignment score threshold.
--pggb-d
Maximum distance for chaining alignments.
--pggb-Q
Prefix used for naming consensus paths.
--pggb-Y
Delimiter used in sequence naming (e.g., "#").
🧠 Notes
All PGGB parameters can be left as default values for general use.
However, it is recommended to adjust them based on genome size, divergence, and ploidy level to achieve optimal results.
🚀 Example Command
sbatch step1-2.sh shuidao \
-p hebhcnormal01 -r shuidao1 -s 1 -e 5 -t 20 \
--pggb-s 5000 --pggb-l 25000 --pggb-p 90 --pggb-c 1 \
--pggb-K 19 --pggb-F 0.001 --pggb-g 30 --pggb-k 23 --pggb-f 0 \
--pggb-B 10M --pggb-n 5 --pggb-j 0 --pggb-e 0 \
--pggb-G 700,900,1100 \
--pggb-P 1,19,39,3,81,1 \
--pggb-O 0.001 --pggb-d 100 \
--pggb-Q Consensus_ --pggb-Y "#"



## Step 2: run PAP

`sbatch --job-name=PAP_pipeline --partition=hebhcnormal01 --nodes=1 --ntasks-per-node=60 --error=%j.err --output=%j.out /public/home/acfurbn1nz/huitian/github/PAP/step3-10.sh --stage1-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --stage2-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --tag ae --threads 60 `


Description:`sbatch --job-name`, `--partition`, `--nodes`, `--ntasks-per-node`, `--error`, and `--output` are parameters used for submitting jobs with `sbatch`. Users should modify these parameters according to the rules of their HPC system.

`--stage1-dir` specifies the input directory for Step 1.

`--stage2-dir` specifies the directory containing all genome `.fa` files that require annotation.

`--tag` refers to the first two letters of the chromosome ID of the reference genome in the VCF file generated in Step 1.

`--threads` specifies the number of threads used by miniprot for transfer annotation.

