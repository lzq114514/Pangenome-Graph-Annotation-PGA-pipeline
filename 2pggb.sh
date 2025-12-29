#!/bin/bash
# ==============================================
# SLURM 参数设置（带默认值）
# ==============================================
#SBATCH --job-name=pggb_and_vg
#SBATCH --partition=hebhcnormal01 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --error=%j.err
#SBATCH --output=%j.out
# 默认参数
DEFAULT_PARTITION="hebhcnormal01"
DEFAULT_REF_PREFIX="putao1"
DEFAULT_WORKDIR="/public/home/acfurbn1nz/huitian/data/gpapipeline/data/putaoworkflow2"
DEFAULT_SINGULARITY_PATH="/public/software/apps/singularity/3.7.3/bin/singularity"
DEFAULT_PGGB_IMAGE="/public/home/acfurbn1nz/pggb-latest.simg"
DEFAULT_PGGB_BIN="/public/home/acfurbn1nz/pggb/pggb"
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
    echo "  -s, --singularity  Singularity path (default: $DEFAULT_SINGULARITY_PATH)"
    echo "  -i, --image        PGGB image path (default: $DEFAULT_PGGB_IMAGE)"
    echo "  -b, --bin          PGGB binary path (default: $DEFAULT_PGGB_BIN)"
    echo "  -t, --threads      Threads per job (default: 20)"
    echo "  -h, --help         Show this help"
    exit 1
}

# 初始化参数
PARTITION="$DEFAULT_PARTITION"
REF_PREFIX="$DEFAULT_REF_PREFIX"
WORKDIR="$DEFAULT_WORKDIR"
SINGULARITY_PATH="$DEFAULT_SINGULARITY_PATH"
PGGB_IMAGE="$DEFAULT_PGGB_IMAGE"
PGGB_BIN="$DEFAULT_PGGB_BIN"
THREADS_PER_JOB=20

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--partition) PARTITION="$2"; shift 2 ;;
        -r|--ref) REF_PREFIX="$2"; shift 2 ;;
        -w|--workdir) WORKDIR="$2"; shift 2 ;;
        -s|--singularity) SINGULARITY_PATH="$2"; shift 2 ;;
        -i|--image) PGGB_IMAGE="$2"; shift 2 ;;
        -b|--bin) PGGB_BIN="$2"; shift 2 ;;
        -t|--threads) THREADS_PER_JOB="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# ==============================================
# 预检查
# ==============================================
cd "$WORKDIR" || { echo "❌ Cannot cd to $WORKDIR"; exit 1; }

[[ -f "$SINGULARITY_PATH" ]] || { echo "❌ Singularity not found: $SINGULARITY_PATH"; exit 1; }
[[ -f "$PGGB_IMAGE" ]] || { echo "❌ PGGB image not found: $PGGB_IMAGE"; exit 1; }
[[ -f "$PGGB_BIN" ]] || { echo "❌ PGGB binary not found: $PGGB_BIN"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not loaded"; exit 1; }
command -v vg >/dev/null || { echo "❌ vg not loaded"; exit 1; }

# PGGB参数
PGGB_PARAMS="-s 5000 -l 25000 -p 90 -c 1 -K 19 -F 0.001 -g 30 -k 23 -f 0 -B 10M -n 12 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -Y \"#\""

# ==============================================
# 步骤1: 索引FASTA文件（串行执行）
# ==============================================
echo -e "\n🔹 Step 1: 索引FASTA文件"
for i in $(seq $START_CHR $END_CHR); do
    fasta_file="${i}all_genomes.fa"
    if [[ -f "$fasta_file" && ! -f "${fasta_file}.fai" ]]; then
        echo "🔄 Indexing: $fasta_file"
        samtools faidx "$fasta_file" || echo "⚠️ Indexing failed (may already exist)"
    fi
done

# ==============================================
# 步骤2: 并行提交PGGB作业（每个染色体独立任务）
# ==============================================
echo -e "\n🔹 Step 2: 提交PGGB作业（并行）"
JOB_IDS=()

for i in $(seq $START_CHR $END_CHR); do
    fasta_file="${i}all_genomes.fa"
    if [[ ! -f "$fasta_file" ]]; then
        echo "⚠️ File not found: $fasta_file"
        continue
    fi

    # 生成每个染色体的独立脚本
    SCRIPT_CONTENT=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=pggb_chr$i
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$THREADS_PER_JOB
#SBATCH --error=chr${i}_pggb.err
#SBATCH --output=chr${i}_pggb.out

echo "🚀 Processing chr$i with $THREADS_PER_JOB threads"
cd "$WORKDIR" || exit 1

# 运行PGGB
$SINGULARITY_PATH exec -B $WORKDIR:/mnt $PGGB_IMAGE $PGGB_BIN \\
    -i "/mnt/$fasta_file" \\
    -o "/mnt/${i}.pggb_out" \\
    $PGGB_PARAMS \\
    --threads $THREADS_PER_JOB

if [[ \$? -ne 0 ]]; then
    echo "❌ PGGB failed for chr$i"
    exit 1
fi

# 运行vg deconstruct
echo "🔄 Running vg deconstruct for chr$i"
GFA_FILE=\$(find "${i}.pggb_out" -name "*.gfa" | head -n 1)
REF_ID=\$(grep -m 1 "^>${REF_PREFIX}#1#chr" "$fasta_file" | cut -d '>' -f 2)

if [[ -z "\$GFA_FILE" || -z "\$REF_ID" ]]; then
    echo "❌ Missing GFA file or reference ID for chr$i"
    exit 1
fi

vg deconstruct -P "\$REF_ID" -a -t $THREADS_PER_JOB -e "\$GFA_FILE" > "${i}.pggb_out/chr${i}.vcf" 2> "${i}.pggb_out/chr${i}.vcf.err"

if [[ \$? -eq 0 ]]; then
    echo "✅ Successfully generated VCF for chr$i"
else
    echo "❌ vg deconstruct failed for chr$i"
    exit 1
fi
EOF
    )

    # 提交作业
    JOB_SCRIPT="pggb_chr${i}.sh"
    echo "$SCRIPT_CONTENT" > "$JOB_SCRIPT"
    JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
    JOB_IDS+=("$JOB_ID")
    echo "📤 Submitted chr$i as job $JOB_ID"
done

# ==============================================
# 输出汇总信息
# ==============================================
echo -e "\n🎉 All jobs submitted!"
echo "➡️ Job IDs: ${JOB_IDS[@]}"
echo "➡️ Monitor with: squeue -u \$USER"
echo "➡️ Output dir: $WORKDIR/*.pggb_out/"
