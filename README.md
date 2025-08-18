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
      <img width="549" height="141" alt="image" src="https://github.com/user-attachments/assets/72292c0b-e1d5-489a-9525-e8c5439ec979" />

 3. Running the pipeline (5 steps):

    Step 1: 
        - Run with:

    `sbatch 1gfchr.sh data/ninanjie`

    - Before running, modify the sbatch parameters in `gffreademapper.py`（bottom of script） and `gfchr.sh`(top of script) to match your job submission system.
    <img width="609" height="234" alt="image" src="https://github.com/user-attachments/assets/0eff2bdf-d286-4cc6-8335-4c53d73f30da" />
    <img width="393" height="204" alt="image" src="https://github.com/user-attachments/assets/5b0a1e8f-b95e-4823-87c6-978b4f5318a0" />



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

    <img width="966" height="204" alt="image" src="https://github.com/user-attachments/assets/2fb9a9d4-be91-4f66-9a30-d271e21bfd6c" />
    <img width="1704" height="33" alt="image" src="https://github.com/user-attachments/assets/53475210-cef9-4d1d-8aeb-05005259224a" />



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

        <img width="1044" height="51" alt="image" src="https://github.com/user-attachments/assets/52a9c32e-d74c-48ce-b929-123977512564" />




    Step 5:
        - Run:

    `python3 5anno.py data/ninanjie`
        
        - Modify the sbatch parameters in the script to match your Linux system, and set the `prefix` parameter to the reference genome prefix followed by "#1#", e.g., for Arabidopsis, set to `ninanjie1#1#`.
<img width="402" height="126" alt="image" src="https://github.com/user-attachments/assets/123a60ae-600a-4465-9d86-23720e3ca1e7" />

        
        - The output will be `anno.gff3`.
