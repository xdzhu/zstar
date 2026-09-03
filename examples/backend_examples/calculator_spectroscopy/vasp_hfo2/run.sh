#!/usr/bin/env bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input="$case_dir/run"
work="$case_dir/work"
dry_run=0
[[ ${1:-} == "--dry-run" ]] && dry_run=1
vasp_command=${VASP_COMMAND:-"mpirun -np 40 vasp_std"}

for file in INCAR INCAR.dfpt KPOINTS POSCAR; do
  [[ -f "$input/$file" ]] || { echo "Missing $input/$file" >&2; exit 3; }
done
missing_potcar=0
if [[ ! -f "$input/POTCAR" ]]; then
  echo "Missing $input/POTCAR (licensed Hf_pv/O PAW data)" >&2
  missing_potcar=1
fi
echo "Case        : $case_dir"
echo "VASP command: $vasp_command"
if [[ $dry_run -eq 1 ]]; then
  [[ $missing_potcar -eq 0 ]] && echo "Dry run passed. No solver was started."
  echo "Stages: cell-relax -> DFPT phonon/BEC/epsilon -> Raman responses"
  exit 0
fi
[[ $missing_potcar -eq 0 ]] || exit 3

command -v zstar >/dev/null 2>&1 || { echo "zstar is not on PATH" >&2; exit 4; }
export OMP_NUM_THREADS=1

relax="$work/00_relax"
mkdir -p "$relax"
cp "$input/INCAR" "$relax/INCAR"
cp "$input/POSCAR" "$relax/POSCAR"
cp "$input/KPOINTS" "$input/POTCAR" "$relax/"
if [[ ! -f "$relax/OUTCAR" ]] || ! grep -q "General timing and accounting" "$relax/OUTCAR"; then
  (cd "$relax" && eval "$vasp_command")
fi

dfpt="$work/01_dfpt"
mkdir -p "$dfpt"
cp "$input/INCAR.dfpt" "$dfpt/INCAR"
cp "$relax/CONTCAR" "$dfpt/POSCAR"
cp "$input/KPOINTS" "$input/POTCAR" "$dfpt/"
if [[ ! -f "$dfpt/OUTCAR" ]] || ! grep -q "General timing and accounting" "$dfpt/OUTCAR"; then
  (cd "$dfpt" && eval "$vasp_command")
fi

spectra="$work/02_spectra"
if [[ ! -f "$spectra/.zstar/spectra.json" ]]; then
  zstar spectra pre --calculator vasp --input-dir "$dfpt" \
    --modes-xml "$dfpt/vasprun.xml" --root "$spectra" --dim 3 --method dfpt
fi
cp "$dfpt/OUTCAR" "$dfpt/vasprun.xml" "$spectra/reference/"
for file in WAVECAR CHGCAR; do
  [[ -f "$dfpt/$file" ]] && cp "$dfpt/$file" "$spectra/reference/$file"
done
zstar spectra run --root "$spectra" --command "$vasp_command" --omp-threads 1
zstar spectra stat --root "$spectra"
zstar spectra post --root "$spectra"

echo "Completed. Generated results are under $work."
