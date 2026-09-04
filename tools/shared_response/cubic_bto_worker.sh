#!/usr/bin/env bash
set -euo pipefail
root=${1:?experiment root}
route=${2:?unified, legacy_bec, or legacy_phonon}
export PATH=/home/zhuxd/Software/anaconda3/envs/icu_copy/bin:$PATH
export PYTHONPATH="$root/code"
export OMP_NUM_THREADS=40 MKL_NUM_THREADS=40 OPENBLAS_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=omp
ulimit -s unlimited
exec python -u "$root/tools/cubic_bto_benchmark.py" run --root "$root" --route "$route"
