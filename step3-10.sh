#!/bin/bash
set -euo pipefail

########################################
# 解析命令行参数
########################################
# Required: --stage1-dir --stage2-dir --tag --threads

while [[ $# -gt 0 ]]; do
    case $1 in
        --stage1-dir) INPUT_DIR_STAGE1="$2"; shift 2 ;;
        --stage2-dir) INPUT_DIR_STAGE2="$2"; shift 2 ;;
        --tag) TAG="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

########################################
# 参数检查
########################################

if [[ -z "${INPUT_DIR_STAGE1:-}" || -z "${INPUT_DIR_STAGE2:-}" || -z "${TAG:-}" || -z "${THREADS:-}" ]]; then
    echo "Usage:"
    echo " pipeline.sh \\"
    echo "   --stage1-dir DIR \\"
    echo "   --stage2-dir DIR \\"
    echo "   --tag TAG \\"
    echo "   --threads THREADS"
    exit 1
fi

########################################
# 自动推导参数
########################################

INPUT_DIR_STAGE1="${INPUT_DIR_STAGE1%/}"
INPUT_DIR_STAGE2="${INPUT_DIR_STAGE2%/}"

PREFIX="$(basename ${INPUT_DIR_STAGE1})1"

SPECIES="$(echo ${PREFIX} | sed 's/[0-9]*$//')"

PROTEIN="${INPUT_DIR_STAGE1}workflow4/miniprotzhushi.fa"

OUTPUT_BASE="${INPUT_DIR_STAGE1}workflow3"

PARTITION="hebhcnormal01"

########################################
# 信息回显
########################################

echo "Stage1 dir: ${INPUT_DIR_STAGE1}"
echo "Stage2 dir: ${INPUT_DIR_STAGE2}"
echo "Prefix: ${PREFIX}"
echo "Tag: ${TAG}"
echo "Species: ${SPECIES}"
echo "Protein: ${PROTEIN}"
echo "Partition: ${PARTITION}"
echo "Threads: ${THREADS}"
echo "4xunz output_base: ${OUTPUT_BASE}"
echo

########################################
# Stage1: 3gfavcf -> patched 4xunz -> 5anno
########################################

echo "===== Stage1: Graph annotation ====="

echo "Step1: 3gfavcf"
python3 3gfavcf.py "${INPUT_DIR_STAGE1}" "${PREFIX}" "${TAG}"
echo "Step1 done."
echo

INPUT_GFF="${INPUT_DIR_STAGE1}/${PREFIX}/${PREFIX}.gff"
INPUT_FNA="${INPUT_DIR_STAGE1}/${PREFIX}/${PREFIX}.fna"

if [[ ! -f "${INPUT_GFF}" ]]; then
    echo "Warning: input_gff not found: ${INPUT_GFF}"
fi

if [[ ! -f "${INPUT_FNA}" ]]; then
    echo "Warning: input_fna not found: ${INPUT_FNA}"
fi

echo "Step2: run 4xunzhaogff (via a patched temporary copy)"

if [[ ! -f "4xunzhaogff.py" ]]; then
    echo "Error: 4xunzhaogff.py not found in $(pwd)"
    exit 1
fi

TMP_SCRIPT="$(mktemp /tmp/4xunzhaogff.XXXXXX.py)"
cp -p 4xunzhaogff.py "${TMP_SCRIPT}"

python3 - <<PYEDIT
import sys,io,re
p = "${TMP_SCRIPT}"
text = open(p, 'r', encoding='utf-8').read()

def repl(var, val, is_path=False, is_expr=False):

    if is_expr:
        new = f"{var} = {val}"
    elif is_path:
        new = f"{var} = Path(\"{val}\")"
    else:
        new = f"{var} = \"{val}\""

    text = re.sub(r"(?m)^" + re.escape(var) + r"\s*=.*$", new, open(p,'r',encoding='utf-8').read(), count=1)
    open(p,'w',encoding='utf-8').write(text)

repl('input_gff', "${INPUT_GFF}")
repl('input_fna', "${INPUT_FNA}")
repl('output_base', "${OUTPUT_BASE}", is_path=True)

fa_expr = 'str(Path(\"${INPUT_DIR_STAGE1}\")) + \"workflow2/*.fa\"'
repl('fa_pattern', fa_expr, is_expr=True)

PYEDIT

python3 "${TMP_SCRIPT}" "${SPECIES}"

RC=$?

rm -f "${TMP_SCRIPT}"

if [[ $RC -ne 0 ]]; then
    echo "4xunzhaogff failed with exit ${RC}"
    exit $RC
fi

echo "Step2 done."
echo

echo "Step3: 5anno"
python3 5anno.py "${INPUT_DIR_STAGE1}"

echo "Step3 done."
echo

echo "===== Stage1 finished ====="
echo

########################################
# Stage2
########################################

echo "===== Stage2: genome annotation ====="

echo "Step4: miniprot"
python3 6miniprot.py \
    --input-dir "${INPUT_DIR_STAGE2}" \
    --protein "${PROTEIN}" \
    --partition "${PARTITION}" \
    --threads "${THREADS}"

echo "Step4 submitted."
echo

########################################
# 等待 gff3 生成
########################################

echo "Waiting for .gff3 files to appear in ${INPUT_DIR_STAGE2} ..."

while true; do
    GFF_COUNT=$(find "${INPUT_DIR_STAGE2}" -type f -name "*.gff3" | wc -l)

    if [[ ${GFF_COUNT} -gt 0 ]]; then
        echo "Detected ${GFF_COUNT} gff3 files."
        break
    fi

    echo "No gff3 yet. Sleeping 60s..."
    sleep 60
done

echo
echo "Step5: sort"
python3 7sort.py "${INPUT_DIR_STAGE2}"

echo "Step5 done."
echo

echo "Step5 done."
echo

echo "Step6: ncbi"
python3 8ncbi.py "${INPUT_DIR_STAGE2}"

echo "Step6 done."
echo

echo "Step7: AS_other"
python3 9AS_other.py "${INPUT_DIR_STAGE2}"

echo "Step7 done."
echo

echo "Step8: pep"
python3 10pep.py "${INPUT_DIR_STAGE2}"

echo "Step8 done."
echo

echo "===== Pipeline finished ====="
