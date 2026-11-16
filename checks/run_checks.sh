#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/run_checks.sh
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Run every check and write checks/RESULTS.md.
#          1. archived run outputs vs the numbers printed in the paper
#          2. rebuilt functions vs the original 2021/2022 functions (optional)
#          3. published figures vs the same figures redrawn from the archive
#          4. rebuilt results (results/) vs the archive and the paper
# Usage:   ADAMANT_DATA=... ADAMANT_ARCHIVE=... ADAMANT_PUBLISHED=... \
#          [ADAMANT_ORIGINAL=... ADAMANT_ORIGINAL_RERUN=...] bash checks/run_checks.sh
#          Run from the repository root after the figure and table scripts.
#          ADAMANT_ORIGINAL (original scripts) enables check 2 (slow, ~1 hour).
# -----------------------------------------------------------------------------
set -euo pipefail
repo=$(pwd)
out=results/checks
mkdir -p "$out"

Rscript checks/check_archive_vs_paper.R | tee "$out/archive_vs_paper.txt"

if [ -n "${ADAMANT_ORIGINAL:-}" ]; then
  Rscript checks/check_equivalence.R | grep -E 'OK|DIFFERENT|equivalence' | tee "$out/equivalence.txt"
fi

# Redraw the figures from the archived run outputs in a scratch copy
redraw=$(mktemp -d)
cp -R R scripts reference "$redraw"/
python3 checks/archive_to_cache.py "$ADAMANT_ARCHIVE" "$redraw/results"
for s in fig2_liver_heat fig3_lambda_hist fig5_eeg_heat fig6_eeg_heat \
         fig7_power_gcv_decreasing fig8_power_gcv_increasing; do
  (cd "$redraw" && Rscript "scripts/$s.R" --cores=1 > /dev/null)
done
python3 checks/check_figures.py "$ADAMANT_PUBLISHED" "$redraw/figures" "$repo" | tee "$out/figures.txt"
rm -rf "$redraw"

python3 checks/check_results.py "$ADAMANT_ARCHIVE" "$repo" ${ADAMANT_ORIGINAL_RERUN:-} | tee "$out/results.txt"
python3 checks/write_results.py "$out" > checks/RESULTS.md
