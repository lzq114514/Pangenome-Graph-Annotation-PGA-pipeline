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
## Step 1: Prepare chromosomes
Run:


`sbatch 1gfchr.sh data/ninanjie`

Description:
Before running, modify the sbatch parameters in gffreademapper.py (bottom of script) and gfchr.sh (top of script) to match your job submission system.
Before running this script, it is recommended to download the official eggNOG databases to the appropriate directory. For example, if you encounter the error
DIAMOND database /public/home/liuzhongqi/miniconda3/envs/zhushi/lib/python3.8/site-packages/data/eggnog_proteins.dmnd not present,
you should create a data directory under
/public/home/liuzhongqi/miniconda3/envs/zhushi/lib/python3.8/site-packages
and then download the required databases from
http://eggnog5.embl.de/download/emapperdb-5.0.2/.
Downloading via command line can be very slow, so it is recommended to download the files locally first and then transfer them to the server.

<img width="300" height="120" alt="image" src="https://github.com/user-attachments/assets/0eff2bdf-d286-4cc6-8335-4c53d73f30da" />

<img width="300" height="120" alt="image" src="https://github.com/user-attachments/assets/5b0a1e8f-b95e-4823-87c6-978b4f5318a0" />

## Step 2: run PAP

`sbatch --job-name=PAP_pipeline --partition=hebhcnormal01 --nodes=1 --ntasks-per-node=60 --error=%j.err --output=%j.out /public/home/acfurbn1nz/huitian/github/PAP/step3-10.sh --stage1-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --stage2-dir /public/home/acfurbn1nz/huitian/public/home/liuzhongqi/pici/twenty --tag ae --threads 60 `


Description:`sbatch --job-name`, `--partition`, `--nodes`, `--ntasks-per-node`, `--error`, and `--output` are parameters used for submitting jobs with `sbatch`. Users should modify these parameters according to the rules of their HPC system.

`--stage1-dir` specifies the input directory for Step 1.

`--stage2-dir` specifies the directory containing all genome `.fa` files that require annotation.

`--tag` refers to the first two letters of the chromosome ID of the reference genome in the VCF file generated in Step 1.

`--threads` specifies the number of threads used by miniprot for transfer annotation.

