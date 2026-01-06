# graph pangenome annotation pipeline (PAP)

"A pggb-based homologous backfill annotation method"

## Table of Contents
1. [Setting up the environment](#setting-up-the-environment)  
2. [Input requirements](#input-requirements)  
3. [Running the pipeline](#running-the-pipeline)  
   - [Step 1: Prepare chromosomes](#step-1-prepare-chromosomes)  
   - [Step 2: Run PGGB](#step-2-run-pggb)  
   - [Step 3: Generate VCF from PGGB](#step-3-generate-vcf-from-pggb)  
   - [Step 4: Unmerge GFF](#step-4-unmerge-gff)  
   - [Step 5: Final annotation](#step-5-final-annotation)
   - [Step 6: Sort annotation files](#step-6-sort-annotation-files)
   - [Step 7: Sort annotation files](#step-7-sort-annotation-files)  

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

<img width="300" height="120" alt="image" src="https://github.com/user-attachments/assets/0eff2bdf-d286-4cc6-8335-4c53d73f30da" />

<img width="300" height="120" alt="image" src="https://github.com/user-attachments/assets/5b0a1e8f-b95e-4823-87c6-978b4f5318a0" />

## Step 2: Run PGGB
Run:


`sbatch 2pggb.sh`

Description:
Edit 2pggb.sh to set the following parameters according to your system:

DEFAULT_PARTITION: job system partition

DEFAULT_REF_PREFIX: prefix of reference genome

DEFAULT_WORKDIR: absolute path to workflow2 (e.g., ninanjieworkflow2)

START_CHR and END_CHR: first and last chromosome numbers of the reference genome

Additional PGGB parameters (e.g., -n, -k, -j) are in the middle of the script.


## Step 3: Generate VCF from PGGB
Run:


`python3 3gfavcf.py data/ninanjie ninanjie1 CP`

Description:

ninanjie1 is the prefix of the reference genome FASTA file used in PGGB

CP is the first two letters of the reference genome chromosome ID (for Arabidopsis chromosome >CP002684.1, use CP)

## Step 4: Unmerge GFF
Run:


`python3 4xunzhaogff.py ninanjie`

Description:
Modify the parameters input_gff and input_fna in the script to point to the reference genome's GFF and FASTA files.The string "ninanjieworkflow3" after output_base should be replaced with "species_workflow3" (i.e., the species name followed by workflow3).

<img width="501" height="147" alt="image" src="https://github.com/user-attachments/assets/d90647d4-bf8c-4edf-a3be-afb5b3a906ea" />


## Step 5: Final annotation
Run:


`python3 5anno.py data/ninanjie`

Description:
Modify sbatch parameters in the script to match your Linux system, and set the prefix parameter to the reference genome prefix followed by #1# (e.g., for Arabidopsis, use ninanjie1#1#).

<img width="300" alt="final annotation prefix" src="https://github.com/user-attachments/assets/123a60ae-600a-4465-9d86-23720e3ca1e7" />

## Step 6: minimap annotation
Run:


`python3 6miniprot.py INPUT_DIR`

Description:INPUT_DIR is the directory containing files to be annotated. It should include all .fa files that need to be annotated.



## Step 7: Sort annotation files
Run:


`python3 7sort.py INPUT_DIR`

Description:
Using 7sort.py, the orphan (gene-unassigned) alternative splicing structures generated in the previous step are organized into a standard annotation format.
