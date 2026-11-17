#!/usr/bin/env bash
# Check and prepare the data root (the external SSD) before any run that
# writes data. It mirrors the checks of riptide::data::requireWritableRoot, so
# a failure here means the simulation would refuse to start as well.
#
# The script never uses sudo, never mounts anything and never creates the data
# root itself: a missing root usually means the SSD is not mounted, and
# creating it would put the data on the system disk. When something is
# missing it prints the exact commands to run by hand.
#
# Usage:
#   scripts/data_root.sh check   [--data-dir DIR] [--need GIB]
#   scripts/data_root.sh prepare [--data-dir DIR] [--need GIB]
#
# Root resolution (same order as resolveDataRoot): --data-dir, then
# $RIPTIDE_DATA_DIR, then /mnt/external_ssd/riptide-sim.
# Environment: RIPTIDE_SSD_MOUNT (default /mnt/external_ssd) is the mount
# point the default root lives on; RIPTIDE_SSD_UUID, if set, is used in the
# suggested mount command.
#
# Exit codes: 0 ok, 2 usage, 3 root missing or SSD not mounted, 4 root not
# writable, 5 root on the system disk, 6 not enough free space.
set -euo pipefail

readonly default_root="/mnt/external_ssd/riptide-sim"
readonly ssd_mount="${RIPTIDE_SSD_MOUNT:-/mnt/external_ssd}"
# Same safety margin as kFreeSpaceMargin in data_paths.hpp.
readonly margin_gib=2

usage() {
    echo "usage: $0 check|prepare [--data-dir DIR] [--need GIB]" >&2
    exit 2
}

print_mount_help() {
    local uuid="${RIPTIDE_SSD_UUID:-<UUID>}"
    cat >&2 <<EOF
To connect the SSD, run by hand (find the UUID with: lsblk -o NAME,SIZE,FSTYPE,UUID):
  sudo mount UUID=${uuid} ${ssd_mount}
  ls -la ${ssd_mount}          # check that it is the right disk
  sudo mkdir -p ${root}        # only the first time
  sudo chown "\$USER":"\$USER" ${root}
EOF
}

[[ $# -ge 1 ]] || usage
command_name="$1"
shift
[[ "${command_name}" == "check" || "${command_name}" == "prepare" ]] || usage

cli_root=""
need_gib=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)
            [[ $# -ge 2 ]] || usage
            cli_root="$2"
            shift 2
            ;;
        --need)
            [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || usage
            need_gib="$2"
            shift 2
            ;;
        *) usage ;;
    esac
done

# An empty value counts as unset, as in resolveDataRoot.
if [[ -n "${cli_root}" ]]; then
    root="${cli_root}"
elif [[ -n "${RIPTIDE_DATA_DIR:-}" ]]; then
    root="${RIPTIDE_DATA_DIR}"
else
    root="${default_root}"
fi
if [[ "${root}" != /* ]]; then
    echo "error: data root must be an absolute path: ${root}" >&2
    exit 2
fi
root="${root%/}"

# 1. The SSD must be mounted when the root lives under its mount point.
if [[ "${root}" == "${ssd_mount}" || "${root}" == "${ssd_mount}"/* ]]; then
    if ! mountpoint -q "${ssd_mount}"; then
        echo "error: ${ssd_mount} is not mounted" >&2
        print_mount_help
        exit 3
    fi
fi

# 2. The root must already exist and be a writable directory.
if [[ ! -d "${root}" ]]; then
    echo "error: data root does not exist: ${root}" >&2
    print_mount_help
    exit 3
fi
if [[ ! -w "${root}" ]]; then
    echo "error: data root is not writable: ${root}" >&2
    echo "  sudo chown \"\$USER\":\"\$USER\" ${root}" >&2
    exit 4
fi

# 3. Never on the device of / (an unmounted SSD leaves an empty mount point there).
if [[ "$(stat -c %d "${root}")" == "$(stat -c %d /)" ]]; then
    echo "error: data root is on the system disk (same device as /): ${root}" >&2
    if [[ "${root}" == "${ssd_mount}"/* ]]; then
        print_mount_help
    else
        echo "  choose a directory on another disk (--data-dir or RIPTIDE_DATA_DIR)" >&2
    fi
    exit 5
fi

# 4. Free space: requested amount plus the safety margin.
free_kib="$(df --output=avail -k "${root}" | tail -n 1 | tr -d ' ')"
required_kib=$(((need_gib + margin_gib) * 1024 * 1024))
if ((free_kib < required_kib)); then
    echo "error: not enough free space on ${root}:" \
        "$((free_kib / 1024 / 1024)) GiB free, ${need_gib} + ${margin_gib} GiB required" >&2
    exit 6
fi

device="$(findmnt -n -o SOURCE --target "${root}")"
echo "data root : ${root}"
echo "device    : ${device} (mounted on $(findmnt -n -o TARGET --target "${root}"))"
echo "free      : $((free_kib / 1024 / 1024)) GiB (required ${need_gib} + ${margin_gib} GiB margin)"

if [[ "${command_name}" == "prepare" ]]; then
    # Subfolders of docs/architecture/data_layout.md; mkdir -p keeps it idempotent.
    mkdir -p "${root}/sim" "${root}/analysis"
    echo "prepared  : ${root}/sim, ${root}/analysis"
fi
echo "ok"
