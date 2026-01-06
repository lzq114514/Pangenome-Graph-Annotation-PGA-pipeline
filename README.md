# Pangenome-Graph-Annotation-PGA-pipeline
"A pggb-based homologous backfill annotation method"
1. Setting up the environment:
    - Use conda to create an environment from pap.yaml.

      `conda env create -f pap.yaml`

2. Input requirements:
    - A folder containing all genome FASTA files (with .fna extension) and corresponding GFF files.
    - For each chromosome ID in the FASTA files, the chromosome must be labeled as "chromosome 1", "chromosome x" for sex chromosomes, etc. (currently only human genome is supported, but the example uses Arabidopsis from the data/ninanjie/fna folder).


      <img width="270" height="70" alt="image" src="https://github.com/user-attachments/assets/72292c0b-e1d5-489a-9525-e8c5439ec979" />

3. Running the pipeline (5 steps):

Step 1: 

- Run with:

`sbatch 1gfchr.sh data/ninanjie`

Description:

Before running, modify the sbatch parameters in `gffreademapper.py`（bottom of script） and `gfchr.sh`(top of script) to match your job submission system.
    

<img width="300" height="120" alt="image" src="https://github.com/user-attachments/assets/0eff2bdf-d286-4cc6-8335-4c53d73f30da" />
    

<img width="200" height="104" alt="image" src="https://github.com/user-attachments/assets/5b0a1e8f-b95e-4823-87c6-978b4f5318a0" />



Step 2:

Description:

- Edit the script `2pggb.sh` to set the following parameters:

DEFAULT_PARTITION: your job system partition

DEFAULT_REF_PREFIX: the prefix of the reference genome

DEFAULT_WORKDIR: the absolute path to workflow2 (e.g., ninanjieworkflow2)

DEFAULT_SINGULARITY_PATH: path to Singularity

DEFAULT_PGGB_IMAGE: path to the downloaded pggb-latest.simg

DEFAULT_PGGB_BIN: path to pggb

START_CHR and END_CHR: the first and last chromosome numbers of the reference genome.

Detailed parameters for PGGB can be found in the middle of this script such as -n -k -j


<img width="470" height="100" alt="image" src="https://github.com/user-attachments/assets/2fb9a9d4-be91-4f66-9a30-d271e21bfd6c" />

    

<img width="1704" height="33" alt="image" src="https://github.com/user-attachments/assets/53475210-cef9-4d1d-8aeb-05005259224a" />



Then run:

`sbatch 2pggb.sh`


Step 3:

Run:

Description:

`python3 3gfavcf.py data/ninanjie ninanjie1 CP`
        
Here, “ninanjie1” is the prefix of the reference genome FASTA file used in pggb, and “CP” is the first two letters of the reference genome chromosome ID (e.g., for Arabidopsis chromosome ">CP002684.1", use "CP").
   

Step 4:
       
Run:

`python3 4xunzhaogff.py ninanjie`

Description:           

Modify the parameters “input_gff” and “input_fna” in the script to point to the reference genome's GFF and FASTA files.

       

<img width="500" height="25" alt="image" src="https://github.com/user-attachments/assets/35487177-4aea-4fe1-ad24-8a86652d6012" />



Step 5:

- Run:
    
`python3 5anno.py data/ninanjie`

Description:
        
Modify the sbatch parameters in the script to match your Linux system, and set the “prefix” parameter to the reference genome prefix followed by "#1#", e.g., for Arabidopsis, set to “ninanjie1#1#”.


<img width="200" height="60" alt="image" src="https://github.com/user-attachments/assets/123a60ae-600a-4465-9d86-23720e3ca1e7" />


Step 6: Miniprot Annotation

Run:

`python miniprot.py --input-dir INPUT_DIR --protein ./workflow4/miniprotzhushi.fa --partition PARTITION --ntasks-per-node NTASKS_PER_NODE --threads THREADS --job-name JOB_NAME --group-workers GROUP_WORKERS`


Description:

--input-dir: Specify the directory containing the genome files to be annotated. Multiple genomes can be included.

--protein: Provide the protein reference database used for annotation. The path should be /workflow4/miniprotzhushi.fa.

Other parameters (e.g., --partition, --threads) can be adjusted according to your system resources.

This step performs homology-based annotation using Miniprot, generating preliminary GFF3 functional annotation files for the input genomes.
        
Step 7: sort Annotation files

Run:

`python 7sort.py INPUT_DIR`

Description:

Using 7sort.py, the orphan (gene-unassigned) alternative splicing structures generated in the previous step were organized into a standard annotation format.


