#!/bin/bash
# The ML scripts here are COPIES of the working versions in the analysis workspace.
# They could not be moved: gin_residual_v2_cv.py is imported by nine other files and
# predict_dataset_gin.py by twelve, all outside this repository.
#
# Duplication is how this project has broken three times already -- predict_dataset_gin
# drifting from the training script (NODE_DIM 32 vs 34, then a stale scaler giving
# uncertainties 1.77x too large), two copies of compute_deltaG.py drifting on atom_ref
# precision, and a pool filter that only one of its two consumers knew about. Each was
# silent. This script cannot prevent the drift, but it makes it loud.
#
#   WORKSPACE=/groups/tluo/FFE_Renzheng ./ml/sync_check.sh
WS="${WORKSPACE:-/groups/tluo/FFE_Renzheng}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
drift=0
for f in active_learning_v3.py gin_residual_v2_cv.py predict_dataset_gin.py \
         tau_calibration.py run_residual.sh; do
    if [ ! -f "$WS/$f" ]; then
        echo "  ?  $f            not found in $WS"
        continue
    fi
    if diff -q "$HERE/$f" "$WS/$f" >/dev/null 2>&1; then
        echo "  ok $f"
    else
        echo "  DRIFT $f  ($(diff "$HERE/$f" "$WS/$f" | grep -c '^[<>]') lines differ)"
        drift=$((drift+1))
    fi
done
# seed_select.py and build_clean_pool.py live here only; the workspace names differ.
for pair in "seed_select.py:seed_select_v2.py" "build_clean_pool.py:build_clean_pool.py"; do
    a="${pair%%:*}"; b="${pair##*:}"
    [ -f "$WS/$b" ] || continue
    if diff -q "$HERE/$a" "$WS/$b" >/dev/null 2>&1; then echo "  ok $a  (workspace: $b)"
    else echo "  DRIFT $a vs workspace $b"; drift=$((drift+1)); fi
done
echo
[ "$drift" -eq 0 ] && echo "in sync" || { echo "$drift file(s) have drifted"; exit 1; }
