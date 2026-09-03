#!/usr/bin/env bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input="$case_dir/run"
work="$case_dir/work"
dry_run=0
[[ ${1:-} == "--dry-run" ]] && dry_run=1

abacus_command=${ABACUS_COMMAND:-abacus}
pyatb_command=${PYATB_COMMAND:-pyatb}
threads=${OMP_NUM_THREADS:-40}
assets=(
  Si_ONCV_PBE-1.0.upf
  C_ONCV_PBE-1.0.upf
  Si_gga_7au_100Ry_2s2p1d.orb
  C_gga_7au_100Ry_2s2p1d.orb
)

for file in STRU INPUT.relax INPUT.scf INPUT.phonon KPT; do
  [[ -f "$input/$file" ]] || { echo "Missing $input/$file" >&2; exit 3; }
done
for file in "${assets[@]}"; do
  [[ -f "$input/assets/$file" ]] || { echo "Missing $input/assets/$file" >&2; exit 3; }
done

echo "Case          : $case_dir"
echo "ABACUS        : $abacus_command"
echo "PYATB         : $pyatb_command"
echo "OpenMP threads: $threads"
if [[ $dry_run -eq 1 ]]; then
  echo "Dry run passed. No solver was started."
  echo "Stages: cell-relax -> BEC -> phonon -> IR/Raman"
  exit 0
fi

command -v zstar >/dev/null 2>&1 || { echo "zstar is not on PATH" >&2; exit 4; }
export OMP_NUM_THREADS="$threads"

relax="$work/00_relax"
mkdir -p "$relax"
cp "$input/STRU" "$relax/STRU"
cp "$input/INPUT.relax" "$relax/INPUT"
cp "$input/KPT" "$relax/KPT"
for file in "${assets[@]}"; do cp "$input/assets/$file" "$relax/$file"; done
if [[ ! -f "$relax/OUT.SIC_RELAX/STRU_ION_D" ]]; then
  (cd "$relax" && eval "$abacus_command")
fi

bec="$work/bec"
mkdir -p "$bec/assets"
cp "$relax/OUT.SIC_RELAX/STRU_ION_D" "$bec/STRU"
cp "$input/INPUT.scf" "$bec/INPUT"
cp "$input/KPT" "$bec/KPT"
for file in "${assets[@]}"; do cp "$input/assets/$file" "$bec/assets/$file"; done
if [[ ! -f "$bec/.zstar/bec.json" ]]; then
  (cd "$bec" && zstar bec pre --stru STRU --pp assets --orb assets \
    --dim 3 --method central --force)
fi
(cd "$bec" && zstar bec run --root . --abacus-command "$abacus_command" \
  --pyatb-command "$pyatb_command" --omp-threads "$threads")
(cd "$bec" && zstar bec stat --root .)
(cd "$bec" && zstar bec post --root .)

phonon="$work/phonon"
mkdir -p "$phonon"
cp "$bec/0.no-move/STRU" "$phonon/STRU"
cp "$input/INPUT.phonon" "$phonon/INPUT"
cp "$input/INPUT.scf" "$phonon/INPUT-scf"
cp "$input/KPT" "$phonon/KPT"
for file in "${assets[@]}"; do cp "$input/assets/$file" "$phonon/$file"; done
if [[ ! -f "$phonon/.zstar/phonon.json" ]]; then
  (cd "$phonon" && zstar phonon pre --root . --stru STRU --dim "1 1 1")
fi
(cd "$phonon" && zstar phonon run --root . --command "$abacus_command" \
  --omp-threads "$threads")
(cd "$phonon" && zstar phonon stat --root .)
(cd "$phonon" && zstar phonon post --root .)
cp "$bec/BORN" "$bec/Z-BORN-symm.out" "$phonon/"

if [[ ! -f "$phonon/raman/.zstar/spectra.json" ]]; then
  copy_args=()
  for file in "${assets[@]}"; do copy_args+=(--copy "$file"); done
  (cd "$phonon" && zstar spectra pre --calculator abacus --kind all \
    --root raman --stru STRU --qpoints qpoints.yaml \
    --born Z-BORN-symm.out --dielectric BORN --modes "4-6" \
    --amplitude 0.02 --copy INPUT-scf --copy KPT "${copy_args[@]}")
fi
(cd "$phonon" && zstar spectra run --root raman \
  --reference "$bec/0.no-move" --abacus-command "$abacus_command" \
  --pyatb-command "$pyatb_command" --omp-threads "$threads")
(cd "$phonon" && zstar spectra stat --root raman)
(cd "$phonon" && zstar spectra post --root raman)

echo "Completed. Generated results are under $work."
