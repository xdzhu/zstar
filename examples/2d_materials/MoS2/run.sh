#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "$ROOT/../../common/run_abacus_case.sh" \
  --case-dir "$ROOT" --work "$ROOT/work" --dim 2 --periodic-axis z "$@"
