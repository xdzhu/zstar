#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
stage=${1:-relax}
case "$stage" in
    relax|fixed_a_control|relax_symmetry_verified) ;;
    *) echo "Choose relax, fixed_a_control, or relax_symmetry_verified" >&2; exit 2 ;;
esac
work="$here/$stage/work"
if [[ -e "$work" ]]; then
    echo "Existing work directory: $work. Preserve it and choose a fresh copy for a new relaxation." >&2
    exit 2
fi
mkdir "$work"
cp -a "$here/$stage/run/." "$work/"
cd "$work"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-40}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-$OMP_NUM_THREADS}
bash -lc "${ABACUS_COMMAND:-abacus}" > run.log 2>&1
