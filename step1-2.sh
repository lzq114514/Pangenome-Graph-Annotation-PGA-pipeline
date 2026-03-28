#!/bin/bash
#SBATCH --job-name=pep_two_stage
#SBATCH --partition=hebhcnormal01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --error=%j.err
#SBATCH --output=%j.out

set -euo pipefail

########################################
# 用法
########################################
usage() {
    echo ""
    echo "Usage:"
    echo "  sbatch $0 <input_dir> [options]"
    echo ""
    echo "Required:"
    echo "  input_dir"
    echo ""
    echo "Pipeline options:"
    echo "  -p, --partition    SLURM partition (default: hebhcnormal01)"
    echo "  -r, --ref          Reference prefix (default: ninanjie1)"
    echo "  -s, --start        Start chromosome index (default: 1)"
    echo "  -e, --end          End chromosome index (default: 5)"
    echo "  -t, --threads      Threads per PGGB job (default: 20)"
    echo ""
    echo "PGGB options:"
    echo "  --pggb-s           (default: 5000)"
    echo "  --pggb-l           (default: 25000)"
    echo "  --pggb-p           (default: 90)"
    echo "  --pggb-c           (default: 1)"
    echo "  --pggb-K           (default: 19)"
    echo "  --pggb-F           (default: 0.001)"
    echo "  --pggb-g           (default: 30)"
    echo "  --pggb-k           (default: 23)"
    echo "  --pggb-f           (default: 0)"
    echo "  --pggb-B           (default: 10M)"
    echo "  --pggb-n           (default: equal to -e, should usually match END_CHR)"
    echo "  --pggb-j           (default: 0)"
    echo "  --pggb-e           (default: 0)"
    echo "  --pggb-G           (default: 700,900,1100)"
    echo "  --pggb-P           (default: 1,19,39,3,81,1)"
    echo "  --pggb-O           (default: 0.001)"
    echo "  --pggb-d           (default: 100)"
    echo "  --pggb-Q           (default: Consensus_)"
    echo "  --pggb-Y           (default: #)"
    echo ""
    echo "Example:"
    echo "  sbatch $0 /data/project -p hebhcnormal01 -r ninanjie1 -s 1 -e 3 -t 40 --pggb-n 3"
    echo ""
    exit 1
}

########################################
# 参数解析
########################################
if [ $# -lt 1 ]; then
    echo "❌ Error: input_dir is required"
    usage
fi

INPUT_DIR="$1"
shift

DEFAULT_PARTITION="hebhcnormal01"
DEFAULT_REF_PREFIX="ninanjie1"

START_CHR=1
END_CHR=5
THREADS_PER_JOB=20

PARTITION="$DEFAULT_PARTITION"
REF_PREFIX="$DEFAULT_REF_PREFIX"

PGGB_S=5000
PGGB_L=25000
PGGB_PCT=90
PGGB_C=1
PGGB_K_CAP=19
PGGB_F=0.001
PGGB_G=30
PGGB_K=23
PGGB_F2=0
PGGB_B="10M"
PGGB_N=""
PGGB_J=0
PGGB_E=0
PGGB_G_LIST="700,900,1100"
PGGB_P_LIST="1,19,39,3,81,1"
PGGB_O=0.001
PGGB_D=100
PGGB_Q="Consensus_"
PGGB_Y="#"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--partition)
            PARTITION="$2"
            shift 2
            ;;
        -r|--ref)
            REF_PREFIX="$2"
            shift 2
            ;;
        -s|--start)
            START_CHR="$2"
            shift 2
            ;;
        -e|--end)
            END_CHR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS_PER_JOB="$2"
            shift 2
            ;;

        --pggb-s)
            PGGB_S="$2"
            shift 2
            ;;
        --pggb-l)
            PGGB_L="$2"
            shift 2
            ;;
        --pggb-p)
            PGGB_PCT="$2"
            shift 2
            ;;
        --pggb-c)
            PGGB_C="$2"
            shift 2
            ;;
        --pggb-K)
            PGGB_K_CAP="$2"
            shift 2
            ;;
        --pggb-F)
            PGGB_F="$2"
            shift 2
            ;;
        --pggb-g)
            PGGB_G="$2"
            shift 2
            ;;
        --pggb-k)
            PGGB_K="$2"
            shift 2
            ;;
        --pggb-f)
            PGGB_F2="$2"
            shift 2
            ;;
        --pggb-B)
            PGGB_B="$2"
            shift 2
            ;;
        --pggb-n)
            PGGB_N="$2"
            shift 2
            ;;
        --pggb-j)
            PGGB_J="$2"
            shift 2
            ;;
        --pggb-e)
            PGGB_E="$2"
            shift 2
            ;;
        --pggb-G)
            PGGB_G_LIST="$2"
            shift 2
            ;;
        --pggb-P)
            PGGB_P_LIST="$2"
            shift 2
            ;;
        --pggb-O)
            PGGB_O="$2"
            shift 2
            ;;
        --pggb-d)
            PGGB_D="$2"
            shift 2
            ;;
        --pggb-Q)
            PGGB_Q="$2"
            shift 2
            ;;
        --pggb-Y)
            PGGB_Y="$2"
            shift 2
            ;;

        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# 默认让 pggb -n 跟随 END_CHR
if [[ -z "$PGGB_N" ]]; then
    PGGB_N="$END_CHR"
fi

########################################
# 自动路径
########################################
OUTPUT_DIR="${INPUT_DIR}/output"
WORKFLOW1_DIR="${INPUT_DIR}workflow1"
WORKDIR="${INPUT_DIR}workflow2"

########################################
# 参数总结
########################################
echo ""
echo "===== PARAMETER SUMMARY ====="
echo "Input dir:        $INPUT_DIR"
echo "Partition:        $PARTITION (default: $DEFAULT_PARTITION)"
echo "Ref prefix:       $REF_PREFIX (default: $DEFAULT_REF_PREFIX)"
echo "Chr range:        $START_CHR - $END_CHR"
echo "Threads/job:      $THREADS_PER_JOB (default: 20)"
echo "Workflow1 dir:    $WORKFLOW1_DIR"
echo "Workflow2 dir:    $WORKDIR"
echo "============================="
echo ""

echo "===== PGGB PARAMETERS ====="
echo "-s $PGGB_S"
echo "-l $PGGB_L"
echo "-p $PGGB_PCT"
echo "-c $PGGB_C"
echo "-K $PGGB_K_CAP"
echo "-F $PGGB_F"
echo "-g $PGGB_G"
echo "-k $PGGB_K"
echo "-f $PGGB_F2"
echo "-B $PGGB_B"
echo "-n $PGGB_N"
echo "-j $PGGB_J"
echo "-e $PGGB_E"
echo "-G $PGGB_G_LIST"
echo "-P $PGGB_P_LIST"
echo "-O $PGGB_O"
echo "-d $PGGB_D"
echo "-Q $PGGB_Q"
echo "-Y $PGGB_Y"
echo "==========================="
echo ""

if [[ "$PGGB_N" != "$END_CHR" ]]; then
    echo "⚠ Warning: PGGB -n ($PGGB_N) is different from END_CHR ($END_CHR)."
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORKFLOW1_DIR"
mkdir -p "$WORKDIR"

echo "-s $PGGB_S -l $PGGB_L -p $PGGB_PCT -c $PGGB_C -K $PGGB_K_CAP -F $PGGB_F -g $PGGB_G -k $PGGB_K -f $PGGB_F2 -B $PGGB_B -n $PGGB_N -j $PGGB_J -e $PGGB_E -G $PGGB_G_LIST -P $PGGB_P_LIST -O $PGGB_O -d $PGGB_D -Q $PGGB_Q -Y $PGGB_Y" > "${WORKDIR}/pggb_params.log"

########################################
# Step 0: 收集 fna 和 gff 到 workflow1
########################################
echo "=== Step 0: collect fna and gff ==="

for sub in "$INPUT_DIR"/*; do
    if [ -d "$sub" ]; then
        name=$(basename "$sub")
        fna="$sub/${name}.fna"
        gff="$sub/${name}.gff"

        [ -f "$fna" ] && cp "$fna" "$WORKFLOW1_DIR/" && echo "✔ $fna" || echo "⚠ missing $fna"
        [ -f "$gff" ] && cp "$gff" "$WORKFLOW1_DIR/" && echo "✔ $gff" || echo "⚠ missing $gff"
    fi
done

########################################
# Step 1: fake emapper
########################################
echo "=== Step 1: fake emapper ==="

for sub in "$INPUT_DIR"/*; do
    if [ -d "$sub" ]; then
        name=$(basename "$sub")
        gff="$sub/${name}.gff"

        if [ -f "$gff" ]; then
            cp "$gff" "$OUTPUT_DIR/${name}_emapper.emapper.gff"
            echo "✔ copied $gff"
        else
            echo "⚠ missing: $gff"
        fi
    fi
done

########################################
# Step 2: 运行 integrated workflow
########################################
echo "=== Step 2: run chrid ==="
python3 chrid.py "$INPUT_DIR"

########################################
# Step 3: 准备并提交 PGGB 作业
########################################
echo "=== Step 3: PGGB ==="

cd "$WORKDIR" || { echo "❌ Cannot cd to $WORKDIR"; exit 1; }

command -v pggb >/dev/null || { echo "❌ pggb not found (conda env not active?)"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }
command -v vg >/dev/null || { echo "❌ vg not found"; exit 1; }

echo "🔹 Index FASTA"
for i in $(seq "$START_CHR" "$END_CHR"); do
    fasta_file="${i}all_genomes.fa"
    if [[ -f "$fasta_file" && ! -f "${fasta_file}.fai" ]]; then
        echo "🔄 Indexing $fasta_file"
        samtools faidx "$fasta_file"
    fi
done

echo "🔹 Submit jobs"
JOB_IDS=()

for i in $(seq "$START_CHR" "$END_CHR"); do
    fasta_file="${i}all_genomes.fa"
    [[ -f "$fasta_file" ]] || continue

    JOB_SCRIPT="pggb_chr${i}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=pggb_chr${i}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${THREADS_PER_JOB}
#SBATCH --error=chr${i}_pggb.err
#SBATCH --output=chr${i}_pggb.out

set -euo pipefail

echo "🚀 Processing chr${i} with ${THREADS_PER_JOB} threads"
cd "$WORKDIR" || exit 1

pggb \\
    -i "$fasta_file" \\
    -o "${i}.pggb_out" \\
    -s "${PGGB_S}" \\
    -l "${PGGB_L}" \\
    -p "${PGGB_PCT}" \\
    -c "${PGGB_C}" \\
    -K "${PGGB_K_CAP}" \\
    -F "${PGGB_F}" \\
    -g "${PGGB_G}" \\
    -k "${PGGB_K}" \\
    -f "${PGGB_F2}" \\
    -B "${PGGB_B}" \\
    -n "${PGGB_N}" \\
    -j "${PGGB_J}" \\
    -e "${PGGB_E}" \\
    -G "${PGGB_G_LIST}" \\
    -P "${PGGB_P_LIST}" \\
    -O "${PGGB_O}" \\
    -d "${PGGB_D}" \\
    -Q "${PGGB_Q}" \\
    -Y "${PGGB_Y}" \\
    --threads "${THREADS_PER_JOB}"

echo "🔄 vg deconstruct chr${i}"

GFA_FILE=\$(find "${i}.pggb_out" -name "*.gfa" | head -n 1)
REF_ID=\$(grep -m 1 "^>${REF_PREFIX}#1#chr" "$fasta_file" | cut -d '>' -f 2)

[[ -n "\$GFA_FILE" && -n "\$REF_ID" ]] || { echo "❌ Missing GFA_FILE or REF_ID"; exit 1; }

vg deconstruct -P "\$REF_ID" -a -t "${THREADS_PER_JOB}" -e "\$GFA_FILE" \\
    > "${i}.pggb_out/chr${i}.vcf" \\
    2> "${i}.pggb_out/chr${i}.vcf.err"

echo "✅ chr${i} done"
EOF

    JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
    JOB_IDS+=("$JOB_ID")
    echo "📤 Submitted chr${i} → job $JOB_ID"
done

echo ""
echo "🎉 Done"
echo "Jobs: ${JOB_IDS[*]}"
echo "Monitor: squeue -u \$USER"
