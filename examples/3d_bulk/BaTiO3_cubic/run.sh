#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export ZSTAR_METHOD=${ZSTAR_METHOD:-forward}
exec "$ROOT/../../common/run_abacus_case.sh" \
  --case-dir "$ROOT" --work "$ROOT/work" --dim 3 --periodic-axis z "$@"
