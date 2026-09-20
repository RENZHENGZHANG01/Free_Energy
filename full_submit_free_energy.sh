#!/bin/bash
# Submit one ORCA free-energy job per molecule.
#
# INPUT_CSV must have a header row and the molecule id (PID) in the FIRST column --
# the same file that scripts/generate_xyz.py consumed, so the ids match the .xyz files
# it wrote. START=2 skips the header; raise it to resume a partially submitted batch.
#
#   INPUT_CSV=data/input_molecules.csv START=2 NUMBER=4000 ./full_submit_free_energy.sh
#
# (The previous default was data/round1.csv with START=29 and NUMBER=1, i.e. a one-off
# resume of a campaign that no longer exists. Those values submitted exactly one job,
# which is not a useful default for anyone picking this up.)

file_path="${INPUT_CSV:-data/input_molecules.csv}"
start="${START:-2}"        # first line to submit; 2 = skip the header
number="${NUMBER:-999999}"  # cap on how many jobs to submit

if [ ! -f "$file_path" ]; then
    echo "ERROR: $file_path not found. Set INPUT_CSV, or create it from a seed set:" >&2
    echo "  python -c \"import pandas as pd; d=pd.read_csv('seeds/seed_v9.csv'); \\" >&2
    echo "             d['PID']='SD'+d['rank'].astype(str); \\" >&2
    echo "             d[['PID','smiles']].to_csv('data/input_molecules.csv',index=False)\"" >&2
    exit 1
fi

echo "Input:  $file_path"
echo "Lines:  from $start, at most $number jobs"

count=0
submitted=0

while IFS=',' read -r MOL _; do
    MOL=$(echo "$MOL" | tr -d '\r')
    count=$((count + 1))

    if [ "$count" -lt "$start" ]; then
        continue
    fi
    [ -z "$MOL" ] && continue

    echo "Submitting free energy job for $MOL"
    qsub -v MOL="$MOL" -N "FE_$MOL" submit_free_energy.csh

    submitted=$((submitted + 1))
    if [ "$submitted" -ge "$number" ]; then
        break
    fi
done < "$file_path"

echo "Submitted $submitted jobs."
