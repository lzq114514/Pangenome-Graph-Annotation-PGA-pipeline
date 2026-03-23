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




## Step 2: run PAP

`sbatch --job-name=PAP_pipeline --partition=hebhcnormal01 --nodes=1 --ntasks-per-node=60 --error=%j.err --output=%j.out /public/home/acfurbn1nz/huitian/github/PAP/step3-10.sh --stage1-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --stage2-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --tag ae --threads 60 `


Description:`sbatch --job-name`, `--partition`, `--nodes`, `--ntasks-per-node`, `--error`, and `--output` are parameters used for submitting jobs with `sbatch`. Users should modify these parameters according to the rules of their HPC system.

`--stage1-dir` specifies the input directory for Step 1.

`--stage2-dir` specifies the directory containing all genome `.fa` files that require annotation.

`--tag` refers to the first two letters of the chromosome ID of the reference genome in the VCF file generated in Step 1.

`--threads` specifies the number of threads used by miniprot for transfer annotation.

