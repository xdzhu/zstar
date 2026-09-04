#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
exec bash "$here/../run_case.sh" "$here" 0 "$@"
