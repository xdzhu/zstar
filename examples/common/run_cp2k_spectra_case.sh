#!/usr/bin/env bash
set -euo pipefail

case_dir=""
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) echo 'Usage: run_cp2k_spectra_case.sh --case-dir DIR [--dry-run]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$case_dir" ]] || { echo 'Missing --case-dir' >&2; exit 2; }
case_dir=$(cd "$case_dir" && pwd)
input="$case_dir/run/input.inp"
root="$case_dir/work"
calculator_command=${CP2K_COMMAND:-cp2k.ssmp}
threads=${OMP_NUM_THREADS:-20}
[[ -f "$input" ]] || { echo "Missing $input" >&2; exit 3; }
echo "Input          : $input"
echo "Work directory : $root"
echo "CP2K command   : $calculator_command"
echo "OpenMP threads : $threads"
if [[ "$dry_run" -eq 1 ]]; then
    cat <<EOF
DRY RUN: no CP2K process will be started.
zstar spectra prepare --calculator cp2k --input "$input" --root "$root" --dim 0
zstar spectra run --root "$root" --command "$calculator_command" --omp-threads "$threads"
zstar spectra stat --root "$root"
zstar spectra collect --root "$root"
EOF
    exit 0
fi
command -v zstar >/dev/null 2>&1 || { echo "zstar is not on PATH." >&2; exit 4; }
data_args=()
if [[ -n "${CP2K_DATA_DIR:-}" ]]; then data_args=(--cp2k-data-dir "$CP2K_DATA_DIR"); fi
zstar spectra prepare --calculator cp2k --input "$input" --root "$root" --dim 0
zstar spectra run --root "$root" --command "$calculator_command" \
    --omp-threads "$threads" "${data_args[@]}"
zstar spectra stat --root "$root"
zstar spectra collect --root "$root"
echo "CP2K IR/Raman spectra completed under $root"
