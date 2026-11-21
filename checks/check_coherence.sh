#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/check_coherence.sh
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Rebuild steps 3 and 4 of the coherence pipeline (coherence/README.md)
#          from the private data with R/coherence.R and with
#          python/adamant/coherence.py, and compare: the files per triple of
#          each with the files the paper used (<data>/04_coherence_SNPs/),
#          and the files of R with those of Python, byte for byte.
#          Not part of `nix flake check`: the data are not in the repository.
#          Step 2 is not checked: the EEG series are not in the data folder.
# Usage:   ADAMANT_DATA=<data folder> bash checks/check_coherence.sh
#          Run from the repository root, inside `nix develop`.
#          The files are rebuilt in a temporary folder that is deleted at
#          the end; only counts and the largest difference are printed.
#          One process at a time, about 4 minutes, 2.7 GB of memory (the
#          genotypes) and 1.5 GB of temporary files.
#          Exit status 1 if a number differs by more than the rounding of
#          the files (15 significant digits), if more than one file differs
#          from the paper's when written correctly rounded, or if R and
#          Python differ.
# -----------------------------------------------------------------------------
set -euo pipefail
data=${ADAMANT_DATA:?set ADAMANT_DATA to the data folder}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

nice Rscript --vanilla checks/check_coherence.R "$data" "$work"
nice python checks/check_coherence.py "$data" "$work"

status=0
if cmp -s <(sort "$work/r/step3.md5") <(sort "$work/python/step3.md5"); then
  echo "R and Python: the $(wc -l < "$work/r/step3.md5" | tr -d ' ') files of step 3 have the same bytes"
else
  echo "R and Python: files of step 3 differ"; status=1
fi
if diff -rq "$work/r/step4" "$work/python/step4" > /dev/null; then
  echo "R and Python: the $(ls "$work/r/step4" | wc -l | tr -d ' ') files of step 4 have the same bytes"
else
  echo "R and Python: files of step 4 differ"; status=1
fi
# Columns of summary.csv: triple, same_bytes, same_bytes_rounded,
# fields_differing, max_difference. The files hold numbers of size up to
# about 10 with 15 significant digits.
awk -F, 'NR > 1 { files++; rounded += $3; if ($5 > worst) worst = $5 }
         END { if (worst > 1e-13 || files - rounded > 1) exit 1 }' "$work/r/summary.csv" \
  || { echo "the rebuilt files differ from the files of the paper"; status=1; }
exit $status
