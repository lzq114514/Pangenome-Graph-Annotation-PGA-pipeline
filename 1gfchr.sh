#!/bin/bash
#SBATCH --job-name=pep
#SBATCH --partition=hebhcnormal01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=60
#SBATCH --error=%j.err
#SBATCH --output=%j.out

set -euo pipefail

#### 参数检查
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

INPUT_DIR="$1"
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ Input directory not found: $INPUT_DIR"
    exit 1
fi

# workflow1 目录（脚本会在这里查找 emapper 注释文件）
WORKFLOW1_DIR="${INPUT_DIR%/}workflow1"

echo "➡ Input: $INPUT_DIR"
echo "➡ Expect workflow1 at: $WORKFLOW1_DIR"

#### 1) 运行 gffreademapper.py 并捕获输出中的 sbatch id（如果有）
echo "🚀 Step 1: Running gffreademapper.py..."
PYOUT=$(python3 gffreademapper.py "$INPUT_DIR" 2>&1 || true)
echo "$PYOUT"

# 提取所有 "Submitted batch job <id>" 中的 id（可能有多个）
readarray -t JOBIDS < <(echo "$PYOUT" | grep -oP 'Submitted batch job \K\d+' || true)

if [ "${#JOBIDS[@]}" -eq 0 ]; then
    echo "⚠️  没有从 gffreademapper.py 输出中捕获到 sbatch job id。"
    echo "   将改为等待 emapper_* 作业和 annotations 文件出现。"
else
    echo "📌 Captured job IDs: ${JOBIDS[*]}"
fi

#### 2) 等待所有捕获到的主 job 和 emapper_* 子作业都完成
echo "⏳ Waiting for captured jobs + emapper_* child jobs to finish..."

WAIT_INTERVAL=60    # 秒，轮询间隔
while true; do
    still_main_running=0

    # 如果有捕获到 JOBIDS，就检查它们
    if [ "${#JOBIDS[@]}" -gt 0 ]; then
        for jid in "${JOBIDS[@]}"; do
            # squeue -j <jid> -h  若有行则还在队列中
            cnt=$(squeue -j "$jid" -h 2>/dev/null | wc -l || true)
            if [ "$cnt" -gt 0 ]; then
                still_main_running=1
                break
            fi
        done
    fi

    # 检查当前用户是否还有以 emapper_ 开头的作业
    # 使用 squeue -u "$USER" -h -o "%j" 然后 grep 计数
    emapper_count=$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -c -E '^emapper_') || emapper_count=0

    # 决策：如果没有 main job 且没有 emapper 作业，则跳出循环
    if [ "$still_main_running" -eq 0 ] && [ "$emapper_count" -eq 0 ]; then
        echo "✅ No main JOBIDs running and no emapper_* jobs running."
        break
    fi

    # 打印状态信息
    if [ "$still_main_running" -eq 1 ]; then
        echo "⏳ Some main job(s) still running... (sleep ${WAIT_INTERVAL}s)"
    fi
    if [ "$emapper_count" -gt 0 ]; then
        echo "⏳ ${emapper_count} emapper_* job(s) still running... (sleep ${WAIT_INTERVAL}s)"
    fi

    sleep "$WAIT_INTERVAL"
done

#### 2.5) 等待 annotations 文件生成（至少一个），超时保护
echo "🔎 Waiting for .emapper.annotations or .annotations files to appear in ${WORKFLOW1_DIR} ..."

if [ ! -d "$WORKFLOW1_DIR" ]; then
    echo "⚠️ workflow1 directory not found: ${WORKFLOW1_DIR}. Creating it (it may be created by upstream script)."
    mkdir -p "$WORKFLOW1_DIR"
fi

ANNOT_TIMEOUT=600   # 秒，最多等 10 分钟（按需调整）
SLEEP_STEP=10
elapsed=0
found=0

while [ "$elapsed" -lt "$ANNOT_TIMEOUT" ]; do
    # 优先查找 *.emapper.annotations，然后再找 *.annotations
    count=$(find "$WORKFLOW1_DIR" -maxdepth 1 -type f \( -name "*.emapper.annotations" -o -name "*.annotations" \) | wc -l || true)
    if [ "$count" -gt 0 ]; then
        echo "✅ Found $count annotation file(s) in ${WORKFLOW1_DIR}."
        found=1
        break
    fi
    echo "⏳ No annotation files yet. waited ${elapsed}s / ${ANNOT_TIMEOUT}s ..."
    sleep "$SLEEP_STEP"
    elapsed=$((elapsed + SLEEP_STEP))
done

if [ "$found" -eq 0 ]; then
    echo "❌ Timeout: no annotation files found in ${WORKFLOW1_DIR} after ${ANNOT_TIMEOUT}s."
    echo "Please check the emapper jobs' outputs and logs in the workflow1 directory."
    exit 1
fi

#### 3) 运行 chrid.py（只有在确认注释文件存在后）
echo "🚀 Step 3: Running chrid.py ..."
python3 chrid.py "$INPUT_DIR"

echo "🎉 chrid.py started."
exit 0
