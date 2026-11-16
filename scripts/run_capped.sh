#!/usr/bin/env bash
# Run a command inside a transient systemd scope with a hard memory cap and no
# swap, so that a runaway simulation or analysis is killed on its own instead
# of driving the whole desktop session out of memory.
#
# Usage: scripts/run_capped.sh [-m MEMORY_MAX] -- command [args...]
#   MEMORY_MAX defaults to 5G (systemd syntax: 512M, 5G, ...).
#
# OMP_NUM_THREADS=1 avoids oversubscription when the command manages its own
# threads (Geant4 MT, TBB).
set -euo pipefail

memory_max="5G"

usage() {
    echo "usage: $0 [-m MEMORY_MAX] -- command [args...]" >&2
    exit 2
}

while getopts ":m:h" opt; do
    case "${opt}" in
        m) memory_max="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))
[[ "${1:-}" == "--" ]] && shift
[[ $# -ge 1 ]] || usage

exec env OMP_NUM_THREADS=1 systemd-run --user --scope \
    -p MemoryMax="${memory_max}" -p MemorySwapMax=0 \
    -- "$@"
