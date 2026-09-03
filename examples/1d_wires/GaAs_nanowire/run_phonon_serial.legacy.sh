#!/usr/bin/env bash
set -euo pipefail

ABACUS_COMMAND=${ABACUS_COMMAND:-"mpirun -np 20 abacus"}
STATE_FILE=${STATE_FILE:-"phonon_stage_state.tsv"}

if [[ ! -f "$STATE_FILE" ]]; then
    printf 'timestamp\tstage\tstate\n' > "$STATE_FILE"
fi

record_state() {
    printf '%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" >> "$STATE_FILE"
}

for stage in disp-*; do
    [[ -d "$stage" ]] || continue
    log="$stage/OUT.PHONON/running_scf.log"
    if [[ -f "$log" ]] && grep -q "Total  Time" "$log"; then
        printf 'SKIP %s\n' "$stage"
        record_state "$stage" SKIPPED_COMPLETE
        continue
    fi
    cp INPUT.phonon "$stage/INPUT"
    cp KPT.phonon "$stage/KPT"
    printf 'START %s\n' "$stage"
    record_state "$stage" RUNNING
    (
        cd "$stage"
        eval "$ABACUS_COMMAND" > abacus.log 2>&1
    )
    if [[ ! -f "$log" ]] || ! grep -q "Total  Time" "$log"; then
        printf 'FAILED %s\n' "$stage" >&2
        record_state "$stage" FAILED
        exit 2
    fi
    printf 'DONE %s\n' "$stage"
    record_state "$stage" COMPLETED
done

printf 'All phonon displacement stages are complete.\n'
