#!/usr/bin/env bash
set -euo pipefail

case_dir=""
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) echo 'Usage: run_cp2k_bec_case.sh --case-dir DIR [--dry-run]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$case_dir" ]] || { echo 'Missing --case-dir' >&2; exit 2; }
case_dir=$(cd "$case_dir" && pwd)
input="$case_dir/run/input.inp"
work="$case_dir/work"
native="$case_dir/native"
cp2k_command=${CP2K_COMMAND:-cp2k.ssmp}
omp_threads=${OMP_NUM_THREADS:-20}
[[ -f "$input" ]] || { echo "Missing $input" >&2; exit 3; }

echo "Case directory : $case_dir"
echo "Input          : $input"
echo "CP2K command   : $cp2k_command"
echo "OpenMP threads : $omp_threads"
if [[ -n "${CP2K_DATA_DIR:-}" ]]; then echo "CP2K data      : $CP2K_DATA_DIR"; fi
if [[ "$dry_run" -eq 1 ]]; then
    cat <<EOF
DRY RUN: no CP2K process will be started.
zstar cp2k-bec prepare --input "$input" --root "$work" --method central --displacement 0.005
zstar cp2k-bec run --root "$work" --cp2k-command "$cp2k_command" --omp-threads "$omp_threads"
zstar cp2k-bec collect --root "$work"
zstar cp2k-bec native --input "$input" --root "$native" --field-strength 1e-4 --cp2k-command "$cp2k_command"
EOF
    exit 0
fi

command -v zstar >/dev/null 2>&1 || { echo "zstar is not on PATH." >&2; exit 4; }
data_args=()
if [[ -n "${CP2K_DATA_DIR:-}" ]]; then data_args=(--data-dir "$CP2K_DATA_DIR"); fi
zstar cp2k-bec prepare --input "$input" --root "$work" \
    --method central --displacement "${ZSTAR_DISPLACEMENT:-0.005}"
zstar cp2k-bec run --root "$work" --cp2k-command "$cp2k_command" \
    --omp-threads "$omp_threads" "${data_args[@]}"
zstar cp2k-bec collect --root "$work"
zstar cp2k-bec native --input "$input" --root "$native" \
    --field-strength "${CP2K_FIELD_STRENGTH:-1e-4}" \
    --cp2k-command "$cp2k_command" --omp-threads "$omp_threads" "${data_args[@]}"
native_apt=$(find "$native" -type f -name '*.data' -print -quit)
[[ -n "$native_apt" ]] || { echo "No native CP2K APT .data file found." >&2; exit 5; }
zstar cp2k-bec compare --zstar-json "$work/cp2k_bec.json" \
    --native-apt "$native_apt" --output "$case_dir/results/comparison.latest.json"
echo "CP2K BEC comparison completed: $case_dir/results/comparison.latest.json"
