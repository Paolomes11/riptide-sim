#!/usr/bin/env bash
# Mount the data SSD on its mount point and prepare the data root, in one
# command. Meant to be run by the user: the mount is the only step that needs
# sudo, and it is the only privileged command here.
#
# The UUID never appears in tracked files. It is read from $RIPTIDE_SSD_UUID,
# then from the untracked file ${XDG_CONFIG_HOME:-~/.config}/riptide-sim/ssd_uuid
# (one line). Mounting by UUID, and refusing anything else on the mount point,
# keeps a different disk from silently becoming the data root.
#
# Running it again when the SSD is already mounted does nothing but the
# final check, so it is safe to call at the start of every session.
#
# Usage: scripts/mount_ssd.sh
# Environment: RIPTIDE_SSD_UUID, RIPTIDE_SSD_MOUNT (default /mnt/external_ssd),
#   RIPTIDE_DATA_DIR (passed on to data_root.sh).
#
# Exit codes: 0 ok, 2 usage or missing/invalid UUID, 3 SSD not connected,
# 7 unsafe state (other disk on the mount point, SSD mounted elsewhere, wrong
# filesystem, mount point missing or not empty), 8 mount failed; any other code
# comes from data_root.sh prepare.
set -euo pipefail

readonly ssd_mount="${RIPTIDE_SSD_MOUNT:-/mnt/external_ssd}"
readonly uuid_file="${XDG_CONFIG_HOME:-${HOME}/.config}/riptide-sim/ssd_uuid"
readonly expected_fstype="ext4"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly script_dir

if [[ $# -ne 0 ]]; then
    echo "usage: $0   (no arguments; see the header for the environment)" >&2
    exit 2
fi

uuid="${RIPTIDE_SSD_UUID:-}"
if [[ -z "${uuid}" && -r "${uuid_file}" ]]; then
    uuid="$(head -n 1 -- "${uuid_file}" | tr -d '[:space:]')"
fi
if [[ -z "${uuid}" ]]; then
    cat >&2 <<EOF
error: SSD UUID unknown. Find it with: lsblk -o NAME,SIZE,FSTYPE,UUID
then either export RIPTIDE_SSD_UUID=<UUID> or save it once:
  mkdir -p "$(dirname -- "${uuid_file}")" && echo <UUID> > "${uuid_file}"
EOF
    exit 2
fi
# Filesystem UUID of ext4: 8-4-4-4-12 hex digits. Anything else is a typo.
if [[ ! "${uuid}" =~ ^[0-9a-fA-F]{8}(-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}$ ]]; then
    echo "error: not a filesystem UUID: ${uuid}" >&2
    exit 2
fi
uuid="${uuid,,}"
readonly device="/dev/disk/by-uuid/${uuid}"

if mountpoint -q "${ssd_mount}"; then
    mounted_source="$(findmnt -n -o SOURCE --mountpoint "${ssd_mount}")"
    mounted_uuid="$(lsblk -n -d -o UUID "${mounted_source}" 2>/dev/null | tr -d '[:space:]')"
    if [[ "${mounted_uuid,,}" != "${uuid}" ]]; then
        echo "error: ${ssd_mount} holds ${mounted_source} (UUID ${mounted_uuid:-?})," \
            "not the data SSD; unmount it first" >&2
        exit 7
    fi
    echo "already mounted: ${mounted_source} on ${ssd_mount}"
else
    if [[ ! -e "${device}" ]]; then
        echo "error: no disk with UUID ${uuid}: is the SSD connected?" >&2
        exit 3
    fi
    fstype="$(lsblk -n -d -o FSTYPE "${device}" | tr -d '[:space:]')"
    if [[ "${fstype}" != "${expected_fstype}" ]]; then
        echo "error: ${device} is '${fstype}', expected ${expected_fstype}" >&2
        exit 7
    fi
    elsewhere="$(findmnt -n -o TARGET --source "UUID=${uuid}" || true)"
    if [[ -n "${elsewhere}" ]]; then
        echo "error: the SSD is already mounted on ${elsewhere}; unmount it or set" \
            "RIPTIDE_SSD_MOUNT" >&2
        exit 7
    fi
    if [[ ! -d "${ssd_mount}" ]]; then
        echo "error: mount point ${ssd_mount} does not exist; create it once with:" >&2
        echo "  sudo mkdir -p ${ssd_mount}" >&2
        exit 7
    fi
    # An unmounted mount point must be empty: files there were written to the
    # system disk while the SSD was missing, and the mount would hide them.
    if [[ -n "$(ls -A -- "${ssd_mount}")" ]]; then
        echo "error: ${ssd_mount} is not mounted but not empty; inspect it first" >&2
        exit 7
    fi
    # sudo can only ask for the password on a terminal; without one (e.g. a
    # tool's captured shell) fail with a clear hint instead of sudo's error.
    if [[ ! -t 0 ]] && ! sudo -n true 2>/dev/null; then
        echo "error: sudo needs a password but there is no terminal to ask it on;" \
            "run this script from a normal terminal" >&2
        exit 8
    fi
    echo "mounting ${device} ($(readlink -f -- "${device}")) on ${ssd_mount} (sudo)"
    if ! sudo mount -t "${expected_fstype}" "UUID=${uuid}" "${ssd_mount}"; then
        echo "error: mount failed" >&2
        exit 8
    fi
    if ! mountpoint -q "${ssd_mount}"; then
        echo "error: ${ssd_mount} is still not a mount point after mount" >&2
        exit 8
    fi
fi

exec "${script_dir}/data_root.sh" prepare
