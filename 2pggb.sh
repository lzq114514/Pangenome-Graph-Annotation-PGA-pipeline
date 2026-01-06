#!/bin/bash
# ==============================================
# SLURM 参数设置
# ==============================================
#SBATCH --job-name=pggb_and_vg
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --error=%j.err
#SBATCH --output=%j.out

# ==============================================
# 默认参数
# ==============================================
DEFAULT_PARTITION="debug"
DEFAULT_REF_PREFIX="putao1"
DEFAULT_WORKDIR="/public/home/yangruo1/zhushi/putaoworkflow2"

START_CHR=1
END_CHR=19

# ==============================================
# 参数解析
# ==============================================
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -p, --partition    SLURM partition (default: $DEFAULT_PARTITION)"
    echo "  -r, --ref          Reference prefix (default: $DEFAULT_REF_PREFIX)"
    echo "  -w, --workdir      Working directory (default: $DEFAULT_WORKDIR)"
    echo "  -t, --threads      Threads per job (default: 20)"
    echo "  -h, --help         Show this help"
    exit 1
}

PARTITION="$DEFAULT_PARTITION"
REF_PREFIX="$DEFAULT_REF_PREFIX"
WORKDIR="$DEFAULT_WORKDIR"
THREADS_PER_JOB=20

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--partition) PARTITION="$2"; shift 2 ;;
        -r|--ref) REF_PREFIX="$2"; shift 2 ;;
        -w|--workdir) WORKDIR="$2"; shift 2 ;;
        -t|--threads) THREADS_PER_JOB="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# ==============================================
# 预检查
# ==============================================
cd "$WORKDIR" || { echo "❌ Cannot cd to $WORKDIR"; exit 1; }

command -v pggb >/dev/null || { echo "❌ pggb not found (conda env not active?)"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }
command -v vg >/dev/null || { echo "❌ vg not found"; exit 1; }

# ==============================================
# PGGB 参数
# ==============================================
PGGB_PARAMS="-s 5000 -l 25000 -p 90 -c 1 -K 19 -F 0.001 -g 30 -k 23 -f 0 \
-B 10M -n 12 -j 0 -e 0 \
-G 700,900,1100 \
-P 1,19,39,3,81,1 \
-O 0.001 -d 100 -Q Consensus_ -Y \"#\""

# ==============================================
# Step 1: 索引 FASTA（串行）
# ==============================================
echo -e "\n🔹 Step 1: Index FASTA"

for i in $(seq $START_CHR $END_CHR); do
    fasta_file="${i}all_genomes.fa"
    if [[ -f "$fasta_file" && ! -f "${fasta_file}.fai" ]]; then
        echo "🔄 Indexing $fasta_file"
        samtools faidx "$fasta_file"
    fi
done

# ==============================================
# Step 2: 提交 PGGB 作业（并行）
# ==============================================
echo -e "\n🔹 Step 2: Submit PGGB jobs"
JOB_IDS=()

for i in $(seq $START_CHR $END_CHR); do
    fasta_file="${i}all_genomes.fa"
    [[ -f "$fasta_file" ]] || continue

    JOB_SCRIPT="pggb_chr${i}.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=pggb_chr$i
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$THREADS_PER_JOB
#SBATCH --error=chr${i}_pggb.err
#SBATCH --output=chr${i}_pggb.out

echo "🚀 Processing chr$i with $THREADS_PER_JOB threads"
cd "$WORKDIR" || exit 1

# -------------------------
# Run PGGB (conda version)
# -------------------------
pggb \\
    -i "$fasta_file" \\
    -o "${i}.pggb_out" \\
    $PGGB_PARAMS \\
    --threads $THREADS_PER_JOB

[[ \$? -eq 0 ]] || { echo "❌ PGGB failed chr$i"; exit 1; }

# -------------------------
# vg deconstruct
# -------------------------
echo "🔄 vg deconstruct chr$i"

GFA_FILE=\$(find "${i}.pggb_out" -name "*.gfa" | head -n 1)
REF_ID=\$(grep -m 1 "^>${REF_PREFIX}#1#chr" "$fasta_file" | cut -d '>' -f 2)

[[ -n "\$GFA_FILE" && -n "\$REF_ID" ]] || exit 1

vg deconstruct -P "\$REF_ID" -a -t $THREADS_PER_JOB -e "\$GFA_FILE" \\
    > "${i}.pggb_out/chr${i}.vcf" \\
    2> "${i}.pggb_out/chr${i}.vcf.err"

echo "✅ chr$i done"
EOF

    JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
    JOB_IDS+=("$JOB_ID")
    echo "📤 Submitted chr$i → job $JOB_ID"
done

# ==============================================
# Summary
# ==============================================
echo -e "\n🎉 All jobs submitted"
echo "Job IDs: ${JOB_IDS[@]}"
echo "Monitor: squeue -u \$USER"


