#!/usr/bin/env bash
set -euo pipefail
root=${1:?experiment root}
case_name=${2:?material case}
export PATH=/home/zhuxd/Software/anaconda3/envs/icu_copy/bin:$PATH
export PYTHONPATH="$root/upload"
export OMP_NUM_THREADS=40 MKL_NUM_THREADS=40 OPENBLAS_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=omp
ulimit -s unlimited
exec python -u "$root/upload/tools/shared_response/complete_eight_systems.py" run --root "$root" --case "$case_name"
