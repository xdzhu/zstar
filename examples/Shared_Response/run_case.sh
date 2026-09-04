#!/usr/bin/env bash
set -euo pipefail

case_dir=${1:?Usage: bash run_case.sh CASE_DIRECTORY DIMENSION}
dimension=${2:?Missing physical dimension}
shift 2
case_dir=$(cd "$case_dir" && pwd)
work=${ZSTAR_WORK:-"$case_dir/work"}
if [[ ! -d "$work" ]]; then
    mkdir -p "$work"
    cp -a "$case_dir/run/." "$work/"
fi
cd "$work"
if [[ ! -f shared_response.json ]]; then
    zstar bec pre --stru STRU --input INPUT --dim "$dimension"
fi
# Paths and MPI/OMP follow zstar config unless command options are supplied.
zstar bec run "$@"
zstar bec stat
for argument in "$@"; do
    if [[ "$argument" == --dry-run ]]; then exit 0; fi
done
zstar bec post
if [[ "$dimension" == 0 ]]; then
    python "$case_dir/../../../tools/shared_response/molecular_validation.py" "$work" \
        --output "$work/molecular_internal_response.json"
else
    zstar dielectric static --dim "$dimension"
fi
