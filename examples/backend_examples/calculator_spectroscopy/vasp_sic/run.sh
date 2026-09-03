#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "$ROOT/../../../common/run_vasp_spectra_case.sh" --case-dir "$ROOT" "$@"
