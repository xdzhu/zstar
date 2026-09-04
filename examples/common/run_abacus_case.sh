#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_abacus_case.sh --case-dir DIR --dim {0,1,2,3} [options]

Options:
  --stage bec|all       Run the BEC/APT workflow, or continue through phonons.
  --work DIR            Work directory (default: CASE_DIR/work).
  --phonon-dim VALUE    Phonopy supercell, e.g. "1 1 2" for the nanowire.
  --periodic-axis AXIS  x, y, or z (default: z).
  --dry-run             Check files and print commands without running a solver.
  -h, --help            Show this help.

Environment:
  ABACUS_COMMAND, PYATB_COMMAND, OMP_NUM_THREADS, ZSTAR_DISPLACEMENT, ZSTAR_METHOD,
  ZSTAR_WORK. Commands may contain launcher arguments, e.g. "mpirun -np 20 abacus".
EOF
}

case_dir=""
dim=""
stage="bec"
work=""
phonon_dim="1 1 1"
periodic_axis="z"
dry_run=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir) case_dir=${2:?missing value for --case-dir}; shift 2 ;;
        --dim) dim=${2:?missing value for --dim}; shift 2 ;;
        --stage) stage=${2:?missing value for --stage}; shift 2 ;;
        --work) work=${2:?missing value for --work}; shift 2 ;;
        --phonon-dim) phonon_dim=${2:?missing value for --phonon-dim}; shift 2 ;;
        --periodic-axis) periodic_axis=${2:?missing value for --periodic-axis}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$case_dir" && -n "$dim" ]] || { usage >&2; exit 2; }
[[ "$stage" == "bec" || "$stage" == "all" ]] || {
    echo "--stage must be bec or all" >&2; exit 2;
}
[[ "$dim" =~ ^[0123]$ ]] || { echo "--dim must be 0, 1, 2, or 3" >&2; exit 2; }

case_dir=$(cd "$case_dir" && pwd)
input_dir="$case_dir/run"
if [[ ! -d "$input_dir" && -f "$case_dir/INPUT" ]]; then
    input_dir="$case_dir"
fi
work=${work:-"${ZSTAR_WORK:-$case_dir/work}"}
if [[ "$work" != /* ]]; then work="$case_dir/$work"; fi
work=$(cd "$(dirname "$work")" && pwd)/$(basename "$work")
abacus_command=${ABACUS_COMMAND:-abacus}
pyatb_command=${PYATB_COMMAND:-pyatb}
omp_threads=${OMP_NUM_THREADS:-20}
displacement=${ZSTAR_DISPLACEMENT:-0.01}
method=${ZSTAR_METHOD:-central}
[[ "$method" == "forward" || "$method" == "central" ]] || {
    echo "ZSTAR_METHOD must be forward or central" >&2; exit 2;
}

for required in INPUT STRU; do
    [[ -f "$input_dir/$required" ]] || {
        echo "Missing $required in $input_dir" >&2; exit 3;
    }
done
[[ -d "$input_dir/assets" ]] || {
    echo "Missing assets/ in $input_dir; the case is not self-contained." >&2; exit 3;
}

echo "Case directory : $case_dir"
echo "Input directory: $input_dir"
echo "Work directory : $work"
echo "Dimension      : $dim"
echo "Stage          : $stage"
echo "ABACUS command : $abacus_command"
echo "PYATB command  : $pyatb_command"
echo "OpenMP threads : $omp_threads"

if [[ "$dry_run" -eq 1 ]]; then
    cat <<EOF
DRY RUN: no files will be generated and no external solver will be started.
1. cp -a "$input_dir/." "$work/"
2. zstar bec pre --stru STRU --input INPUT --pp assets --orb assets --dim $dim \\
     --method $method --displacement $displacement --force
3. zstar workflow run --root . --dimensionality $dim --omp-threads $omp_threads \\
     --abacus-command "$abacus_command" --pyatb-command "$pyatb_command"
4. zstar bec post --root .
5. Reuse shared Gamma forces for --phonon-dim "1 1 1"; otherwise prepare
   a separate phonon/ directory, using INPUT.phonon when provided.
EOF
    exit 0
fi

command -v zstar >/dev/null 2>&1 || {
    echo "zstar is not on PATH; install the package or activate its environment." >&2
    exit 4
}

if [[ ! -f "$work/.zstar-example-seeded" ]]; then
    mkdir -p "$work"
    cp -a "$input_dir/." "$work/"
    printf 'Seeded from %s\n' "$input_dir" > "$work/.zstar-example-seeded"
fi

cd "$work"
if [[ ! -f .zstar-bec-prepared ]]; then
    zstar bec pre --stru STRU --input INPUT --pp assets --orb assets --dim "$dim" \\
        --method "$method" --displacement "$displacement" --force
    touch .zstar-bec-prepared
fi

zstar workflow run --root . --dimensionality "$dim" --omp-threads "$omp_threads" \\
    --abacus-command "$abacus_command" --pyatb-command "$pyatb_command"
zstar workflow status --root .
zstar bec post --root .

if [[ "$stage" == "all" ]]; then
    if [[ -f shared_response.json && "$phonon_dim" == "1 1 1" ]]; then
        zstar phonon post --stru STRU --physical-dim "$dim"
        echo "Reused the shared BEC/Gamma response ensemble."
    else
        # Keep supercell force files out of a shared Gamma ensemble.
        if [[ -f shared_response.json ]]; then
            mkdir -p phonon
            if [[ ! -f phonon/phonopy_disp.yaml ]]; then
                cp STRU phonon/STRU
                cp -a assets phonon/
                [[ ! -f KPT ]] || cp KPT phonon/KPT
                cp INPUT phonon/INPUT
            fi
            cd phonon
        fi
        if [[ -f "$input_dir/INPUT.phonon" ]]; then
            cp "$input_dir/INPUT.phonon" INPUT
        fi
        if [[ ! -f phonopy_disp.yaml ]]; then
            zstar ph --stru STRU --dim "$phonon_dim"
        fi
        zstar phonon run --root . --command "$abacus_command" --omp-threads "$omp_threads"
        zstar postph --stru STRU --physical-dim "$dim"
        zstar phonon stat --root . || true
    fi
fi

echo "ZStar example stage '$stage' completed or resumed successfully."
