#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
work=${ZSTAR_WORK:-"$here/work"}
if [[ ! -d "$work" ]]; then
    mkdir -p "$work"
    cp -a "$here/run/." "$work/"
fi
cd "$work"
if [[ ! -d 0.no-move ]]; then
    zstar bec pre --stru STRU --input INPUT --dim 2 --ensemble cartesian \
        --method forward --displacement 0.01
fi
zstar bec run "$@"
zstar bec stat
for argument in "$@"; do
    if [[ "$argument" == --dry-run ]]; then exit 0; fi
done
zstar bec post --dim 2 --method forward
