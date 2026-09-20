#!/bin/bash
# Recompute the size-independent ΔG residual for the ML workspace's dataset.
#
# compute_residual_deltaG.py now lives with the ORCA pipeline (Renzheng/scripts/,
# where it is Step 6) so there is only ONE copy of the FS5 baseline. It defaults to
# the dataset next to itself, so point it at this directory's data explicitly.
#
# Fast (Ridge on ~1.3k rows) — no need for qsub, just run it on the login node.
#
#   ./run_residual.sh                       # -> final_data_with_residual_deltaG.csv
#   IN=final_data_with_deltaG_1.csv OUT=final_data_with_residual_deltaG_1.csv ./run_residual.sh
set -e
export PATH="/users/rzhang4/.conda/envs/rdkit_env/bin:$PATH"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IN="${IN:-final_data_with_deltaG.csv}"
OUT="${OUT:-final_data_with_residual_deltaG.csv}"

export RESIDUAL_INPUT="$DIR/$IN"
export RESIDUAL_OUTPUT="$DIR/$OUT"
export RESIDUAL_PLOTDIR="$DIR/residual_plots"

echo "  in  : $RESIDUAL_INPUT"
echo "  out : $RESIDUAL_OUTPUT"
echo "  plots: $RESIDUAL_PLOTDIR"
echo

python -u "$DIR/Renzheng/scripts/compute_residual_deltaG.py"
