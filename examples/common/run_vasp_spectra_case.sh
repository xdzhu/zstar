#!/usr/bin/env bash
set -euo pipefail

case_dir=""
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) echo 'Usage: run_vasp_spectra_case.sh --case-dir DIR [--dry-run]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$case_dir" ]] || { echo 'Missing --case-dir' >&2; exit 2; }
case_dir=$(cd "$case_dir" && pwd)
input="$case_dir/run"
root="$case_dir/work"
calculator_command=${VASP_COMMAND:-vasp_std}
threads=${OMP_NUM_THREADS:-20}
missing=()
for file in INCAR POSCAR KPOINTS POTCAR vasprun.xml; do
    [[ -f "$input/$file" ]] || missing+=("$file")
done
echo "Input directory : $input"
echo "Work directory  : $root"
echo "VASP command    : $calculator_command"
if [[ "$dry_run" -eq 1 ]]; then
    [[ ${#missing[@]} -eq 0 ]] && echo "All VASP inputs are present." || echo "Missing licensed/local inputs: ${missing[*]}"
    echo "No VASP process will be started."
    exit 0
fi
if [[ ${#missing[@]} -ne 0 ]]; then
    echo "Missing VASP inputs in $input: ${missing[*]}" >&2
    echo "POTCAR must be obtained under the user's VASP license." >&2
    exit 3
fi
command -v zstar >/dev/null 2>&1 || { echo "zstar is not on PATH." >&2; exit 4; }
zstar spectra prepare --calculator vasp --input-dir "$input" \
    --modes-xml "$input/vasprun.xml" --root "$root" --dim 3 --method dfpt
zstar spectra run --root "$root" --command "$calculator_command" --omp-threads "$threads"
zstar spectra stat --root "$root"
zstar spectra collect --root "$root"
echo "VASP IR/Raman spectra completed under $root"
