#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
dry_run=0
for arg in "$@"; do
  [[ "$arg" == "--dry-run" ]] && dry_run=1
done

"$ROOT/../../common/run_abacus_case.sh" \
  --case-dir "$ROOT" --work "$ROOT/work" --dim 2 --stage all \
  --phonon-dim "4 4 1" --periodic-axis z "$@"

if [[ "$dry_run" -eq 1 ]]; then
  cat <<'EOF'
6. zstar dielectric static --qpoints qpoints.yaml \
     --born Z-BORN-symm.out --dielectric BORN --dim 2
7. zstar dielectric freq --qpoints qpoints.yaml \
     --born Z-BORN-symm.out --dielectric BORN --dim 2 \
     --broadening 8 --max-frequency 1600
EOF
  exit 0
fi

cd "$ROOT/work"
zstar dielectric static --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2
zstar dielectric freq --qpoints qpoints.yaml \
  --born Z-BORN-symm.out --dielectric BORN --dim 2 \
  --broadening 8 --max-frequency 1600
