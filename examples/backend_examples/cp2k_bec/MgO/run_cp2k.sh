#!/usr/bin/env bash
set -euo pipefail

run_root="${1:?usage: run_cp2k.sh RUN_ROOT}"
python_bin="${ZSTAR_PYTHON:-python}"
cp2k_command="${CP2K_COMMAND:-cp2k.ssmp}"
omp_threads="${OMP_NUM_THREADS:-20}"
data_args=()

if [[ -n "${CP2K_DATA_DIR:-}" ]]; then
  data_args=(--data-dir "${CP2K_DATA_DIR}")
fi

cd "${run_root}"
"${python_bin}" -m zstar.cli cp2k-bec run \
  --root work \
  --cp2k-command "${cp2k_command}" \
  --omp-threads "${omp_threads}" \
  "${data_args[@]}"
