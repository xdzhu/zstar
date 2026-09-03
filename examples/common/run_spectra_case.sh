#!/usr/bin/env bash
set -euo pipefail

case_dir=""
dim=""
work=""
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --dim) dim=${2:?missing value for --dim}; shift 2 ;;
        --work) work=${2:?missing value for --work}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) echo 'Usage: run_spectra_case.sh --case-dir DIR --dim {0,1,2,3} [--work DIR] [--dry-run]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$case_dir" && -n "$dim" ]] || { echo 'Missing --case-dir or --dim' >&2; exit 2; }
case_dir=$(cd "$case_dir" && pwd)
work=${work:-"$case_dir/work"}
if [[ "$work" != /* ]]; then work="$case_dir/$work"; fi
abacus_command=${ABACUS_COMMAND:-abacus}
pyatb_command=${PYATB_COMMAND:-pyatb}
omp_threads=${OMP_NUM_THREADS:-20}

if [[ "$dry_run" -eq 1 ]]; then
    "$case_dir/../../common/run_abacus_case.sh" --case-dir "$case_dir" \
        --work "$work" --dim "$dim" --stage all --dry-run
    cat <<EOF
6. zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
     --dielectric BORN --dim $dim --outdir ir_spectrum
7. zstar raman prepare --stru STRU --qpoints qpoints.yaml --outdir raman
8. zstar raman run --raman-dir raman --reference 0.no-move \
     --qpoints qpoints.yaml --dim $dim
DRY RUN: no spectroscopy solver or post-processing command was started.
EOF
    exit 0
fi

"$case_dir/../../common/run_abacus_case.sh" --case-dir "$case_dir" \
    --work "$work" --dim "$dim" --stage all
cd "$work"
zstar ir --qpoints qpoints.yaml --born Z-BORN-symm.out \
    --dielectric BORN --dim "$dim" --outdir ir_spectrum
zstar raman prepare --stru STRU --qpoints qpoints.yaml --outdir raman
zstar raman run --raman-dir raman --reference 0.no-move \
    --qpoints qpoints.yaml --dim "$dim" \
    --abacus-command "$abacus_command" --pyatb-command "$pyatb_command"
echo "IR/Raman spectra completed under $work"
