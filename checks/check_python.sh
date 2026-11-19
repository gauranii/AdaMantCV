#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/check_python.sh
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Build every table and figure with the R scripts (scripts/) and
#          with the Python scripts (python/scripts/), time both, and compare
#          the two by code (checks/check_python.py): every cell of the
#          tables and the numbers every figure plots.
# Usage:   bash checks/check_python.sh archive|rerun <work dir> [PYTHON.md]
#          [folder of the comparison tables]
#          Run from the repository root, inside `nix develop`.
#          archive: both from the archived runs (the defaults of the
#                   scripts); this is what `nix flake check` runs.
#          rerun:   both recomputed from the data (needs ADAMANT_DATA for the
#                   EEG outputs): Table 4 and Figures 1, 2, 3, 7, 8 in full,
#                   Tables 1 to 3 with 2 replications per setting and
#                   Figures 5 and 6 with 3 replications per triple.
#          The outputs go to <work dir>/<mode>/{r,python}; finished outputs
#          are kept, so an interrupted run continues. The comparison is
#          written to python_<mode>.csv in the folder given last (checks/ by
#          default), and both modes together to the Markdown file given
#          (checks/PYTHON.md by default).
#          CORES (default 4) is passed to the scripts as --cores; IMPLS
#          (default "r python") chooses the implementations to run.
#          A run stops when its processes use more than MAX_RSS_MB (3072).
#          Exit status 1 if a comparison fails.
# -----------------------------------------------------------------------------
set -euo pipefail
mode=$1
work=$2
report=${3:-checks/PYTHON.md}
tables=${4:-checks}  # folder of python_<mode>.csv
cores=${CORES:-4}
impls=${IMPLS:-r python}
max_rss_mb=${MAX_RSS_MB:-3072}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

outputs="table1_type1_error table2_power_fixed table3_power_random table4_eeg_pvalues
fig1_liver_cv fig2_liver_heat fig3_lambda_hist fig5_eeg_heat fig6_eeg_heat
fig7_power_gcv_decreasing fig8_power_gcv_increasing"

# Options of one script in the given mode.
options() {
  [ "$mode" = rerun ] || return 0
  case $1 in
    fig1_*) ;;
    table[123]_*) echo "--rerun --runs=2 --cores=$cores" ;;
    fig[56]_*) echo "--rerun --runs=3 --cores=$cores" ;;
    *) echo "--rerun --cores=$cores" ;;
  esac
}

# Resident memory in MB of a process and its descendants.
tree_rss_mb() {
  ps -axo pid=,ppid=,rss= 2>/dev/null | awk -v root="$1" '
    { parent[$1] = $2; rss[$1] = $3 }
    END { for (p in parent) { q = p; while (q != root && q in parent && q != 1) q = parent[q]
                              if (q == root) total += rss[p] }
          printf "%d", total / 1024 }'
}

# Run a command, record its time, stop it above the memory limit.
timed() {
  local impl=$1 name=$2; shift 2
  local dir=$work/$mode/$impl
  mkdir -p "$dir/logs"
  [ -e "$dir/logs/$name.done" ] && return 0
  local start=$EPOCHREALTIME
  ADAMANT_OUTPUT=$dir ADAMANT_PLOT_DATA=$dir/plot_data "$@" > "$dir/logs/$name.log" 2>&1 &
  local pid=$! peak=0 rss
  while kill -0 $pid 2>/dev/null; do
    rss=$(tree_rss_mb $pid || true); rss=${rss:-0}
    [ "$rss" -gt "$peak" ] && peak=$rss
    if [ "$rss" -gt "$max_rss_mb" ]; then
      pkill -P $pid || true; kill $pid || true
      echo "$impl $name: stopped at $rss MB" >&2; exit 2
    fi
    sleep 1
  done
  wait $pid || { echo "$impl $name failed, see $dir/logs/$name.log" >&2; exit 1; }
  awk -v a="$start" -v b="$EPOCHREALTIME" -v p="$peak" 'BEGIN { printf "%.1f,%d\n", b - a, p }' \
    > "$dir/logs/$name.done"
  echo "$impl $name: $(cat "$dir/logs/$name.done") (s, peak MB)"
}

for name in $outputs; do
  for impl in $impls; do
    if [ "$impl" = r ]; then
      timed r "$name" Rscript checks/plot_data.R "$name" $(options "$name")
    else
      timed python "$name" python python/scripts/"$name".py $(options "$name")
    fi
  done
done

[ "$impls" = "r python" ] || exit 0
status=0
python checks/check_python.py "$mode" "$work/$mode/r" "$work/$mode/python" \
  "$tables/python_$mode.csv" || status=1
other=$([ "$mode" = archive ] && echo rerun || echo archive)
[ -e "$tables/python_$other.csv" ] || cp "checks/python_$other.csv" "$tables/"
python checks/write_python.py "$tables/python_archive.csv" "$tables/python_rerun.csv" > "$report"
exit $status
