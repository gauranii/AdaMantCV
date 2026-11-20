#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/versions.sh
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Print the versions of the software the flake provides: the nixpkgs
#          revision of flake.lock, R and its packages, Python and the packages
#          of pyEnv in flake.nix, and the BLAS and LAPACK library.
# Usage:   bash checks/versions.sh            print the versions
#          bash checks/versions.sh --check    compare them with the section
#                                             "Software versions" of README.md
#          Run from the repository root, inside `nix develop`.
# Outputs: prints the versions; with --check exit status 1 if README.md differs.
# -----------------------------------------------------------------------------
set -euo pipefail

versions() {
  python - <<'PY'
import json
from importlib.metadata import version
import platform

locked = json.load(open("flake.lock"))["nodes"]["nixpkgs"]["locked"]
print(f"nixpkgs  {locked['owner']}/{locked['repo']} {locked['rev']}")
print()
print(f"Python   {platform.python_version()}")
for name in ["numpy", "scipy", "pandas", "matplotlib", "numba"]:  # pyEnv of flake.nix
    print(f"  {name:12} {version(name)}")
PY
  echo
  Rscript -e '
    cat(sprintf("R        %s.%s\n", R.version$major, R.version$minor))
    installed <- installed.packages()
    installed <- installed[is.na(installed[, "Priority"]) | installed[, "Priority"] != "base", ]
    installed <- installed[order(tolower(rownames(installed))), ]
    cat(sprintf("  %-12s %s\n", rownames(installed), installed[, "Version"]), sep = "")'
  echo
  echo "OpenBLAS $(echo "$ADAMANT_BLAS" | sed -E 's|.*-openblas-([0-9.]+)/.*|\1|')"
}

if [ "${1:-}" = --check ]; then
  diff <(versions) <(awk '/^## Software versions/ { section = 1 }
                          section && /^```/ { block++; next }
                          section && block == 1 { print }' README.md)
else
  versions
fi
