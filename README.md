# Pangenome-Graph-Annotation-PGA-pipeline
"A pggb-based homologous backfill annotation method"
1. Setting up the environment:
    - Use conda to create an environment from pga.yaml.

      `conda env create -f pga.yaml`
    - Use Singularity to download the pggb tool

      `singularity pull docker://ghcr.io/pangenome/pggb:latest`
    - Clone the pggb repository

      `git clone --recursive https://github.com/pangenome/pggb.git`
 2. Input requirements:
    - A folder containing all genome FASTA files (with .fna extension) and corresponding GFF files.
    - For each chromosome ID in the FASTA files, the chromosome must be labeled as "chromosome 1", "chromosome x" for sex chromosomes, etc. (currently only human genome is supported, but the example uses Arabidopsis from the data/ninanjie/fna folder).
 3. Running the pipeline (5 steps):

    Step 1: 
        - Run with:

    `sbatch 1gfchr.sh data/ninanjie`
        - Before running, modify the sbatch parameters in `gffreademapper.py`（bottom of script） and `gfchr.sh`(top of script) to match your job submission system.


    Step 2:
        - Edit the script `2pggb.sh` to set the following parameters:

    DEFAULT_PARTITION: your job system partition

    DEFAULT_REF_PREFIX: the prefix of the reference genome

    DEFAULT_WORKDIR: the absolute path to workflow2 (e.g., ninanjieworkflow2)

    DEFAULT_SINGULARITY_PATH: path to Singularity

    DEFAULT_PGGB_IMAGE: path to the downloaded pggb-latest.simg

    DEFAULT_PGGB_BIN: path to pggb

    START_CHR and END_CHR: the first and last chromosome numbers of the reference genome.

    Detailed parameters for PGGB can be found in the middle of this script such as -n -k -j
        - Then run:

    `sbatch 2pggb.sh`


    Step 3:
        - Run:
   `python3 3gfavcf.py data/ninanjie ninanjie1 CP`
        - Here, `ninanjie1` is the prefix of the reference genome FASTA file used in pggb, and `CP` is the first two letters of the reference genome chromosome ID (e.g., for Arabidopsis chromosome ">CP002684.1", use "CP").
   

    Step 4:
        - Run:

    `python3 4xunzhaogff.py ninanjie`
        - Modify the parameters `input_gff` and `input_fna` in the script to point to the reference genome's GFF and FASTA files.


    Step 5:
        - Run:

    `python3 5anno.py data/ninanjie`
        - Modify the sbatch parameters in the script to match your Linux system, and set the `prefix` parameter to the reference genome prefix followed by "#1#", e.g., for Arabidopsis, set to `ninanjie1#1#`.
        - The output will be `anno.gff3`.
