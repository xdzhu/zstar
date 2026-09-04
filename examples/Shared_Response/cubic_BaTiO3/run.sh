#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
work=${ZSTAR_WORK:-"$here/work"}
if [[ ! -d "$work" ]]; then
    mkdir -p "$work"
    cp -a "$here/run/." "$work/"
fi
cd "$work"
if [[ ! -f shared_response.json ]]; then
    zstar bec pre --stru STRU -i INPUT --pp assets --orb assets
fi
zstar bec run --mp-density 0.08 "$@"
zstar bec stat
for argument in "$@"; do
    if [[ "$argument" == --dry-run ]]; then exit 0; fi
done
zstar bec post
zstar phonon irrep
printf '%s\n' 'BEC and Gamma modes collected. Cubic BTO has a soft mode; no stable static phonon dielectric constant is reported.'
