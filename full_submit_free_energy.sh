#!/bin/bash
file_path="data/PI1070.csv"   # CSV 第一列是 molecule ID
start=781                 # 从第几行开始（如跳过表头就设为 2）
number=100         # 最多提交多少个

count=0
submitted=0

while IFS=',' read -r MOL _; do
    MOL=$(echo "$MOL" | tr -d '\r')
    count=$((count + 1))

    if [ "$count" -lt "$start" ]; then
        continue
    fi

    echo "Submitting free energy job for $MOL"
    qsub -v MOL="$MOL" -N "FE_$MOL" submit_free_energy.csh

    submitted=$((submitted + 1))
    if [ "$submitted" -ge "$number" ]; then
        break
    fi
done < "$file_path"
