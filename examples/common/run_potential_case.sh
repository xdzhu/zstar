#!/usr/bin/env bash
set -euo pipefail

case_dir=""
mode="slab"
cube=""
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --mode) mode=${2:?missing value for --mode}; shift 2 ;;
        --cube) cube=${2:?missing value for --cube}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) echo 'Usage: run_potential_case.sh --case-dir DIR --mode mos2|in2se3|directional [--cube FILE] [--dry-run]'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$case_dir" ]] || { echo 'Missing --case-dir' >&2; exit 2; }
case_dir=$(cd "$case_dir" && pwd)
input_dir="$case_dir/run"
cube=${cube:-"${ZSTAR_CUBE:-$case_dir/work/ElecStaticPot.cube}"}
case "$mode" in
    mos2)
        [[ -f "$input_dir/INPUT" && -f "$input_dir/STRU" ]] || { echo "Missing INPUT/STRU in $input_dir" >&2; exit 3; }
        asset_root="$input_dir/assets"
        [[ -d "$asset_root" ]] || asset_root="$input_dir"
        compgen -G "$asset_root/*.upf" >/dev/null || { echo "Missing pseudopotential files in $asset_root" >&2; exit 3; }
        compgen -G "$asset_root/*.orb" >/dev/null || { echo "Missing orbital files in $asset_root" >&2; exit 3; }
        extra=(--center-slab) ;;
    in2se3)
        [[ -f "$input_dir/INPUT" && -f "$input_dir/STRU" ]] || { echo "Missing INPUT/STRU in $input_dir" >&2; exit 3; }
        asset_root="$input_dir/assets"
        [[ -d "$asset_root" ]] || asset_root="$input_dir"
        compgen -G "$asset_root/*.upf" >/dev/null || { echo "Missing pseudopotential files in $asset_root" >&2; exit 3; }
        compgen -G "$asset_root/*.orb" >/dev/null || { echo "Missing orbital files in $asset_root" >&2; exit 3; }
        extra=(--polar-arrow auto) ;;
    directional)
        extra=(--direction a+b --direction a-b --direction-bins 160 --direction-method linear --direction-samples 72 72 --direction-smooth 0.15 --mirror-test) ;;
    *) echo '--mode must be mos2, in2se3, or directional' >&2; exit 2 ;;
esac
echo "Case directory : $case_dir"
echo "Input directory: $input_dir"
echo "Cube input     : $cube"
if [[ "$dry_run" -eq 1 ]]; then
    cat <<EOF
DRY RUN: no files will be generated and no solver will be started.
1. Provide an electrostatic-potential cube from a converged ABACUS calculation.
2. zstar pot --cube "$cube" --axes z --plane xy --plane-average --tile 5 5 \
     --vacuum-level --vacuum-sides --vacuum-window 0.75 \
     ${extra[*]} --outdir "$case_dir/work/potential"
EOF
    exit 0
fi
command -v zstar >/dev/null 2>&1 || { echo 'zstar is not on PATH.' >&2; exit 4; }
[[ -f "$cube" ]] || {
    echo "Missing electrostatic-potential cube: $cube" >&2
    echo 'Run the DFT SCF first or pass --cube PATH.' >&2
    exit 5
}
mkdir -p "$case_dir/work/potential"
zstar pot --cube "$cube" --axes z --plane xy --plane-average --tile 5 5 \
    --vacuum-level --vacuum-sides --vacuum-window 0.75 \
    "${extra[@]}" --outdir "$case_dir/work/potential"
echo "Potential analysis completed under $case_dir/work/potential"
