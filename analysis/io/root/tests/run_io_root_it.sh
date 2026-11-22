#!/usr/bin/env bash
# Integration test of riptide::io_root on the 1- and 2-thread runs of the
# riptide_sim integration test (ctest fixture sim_integration). The reader
# only reads the runs; its one write is a small file with a column removed,
# under a reserved analysis directory that is recreated every time.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the integration runs.
#
# Usage: run_io_root_it.sh IO_ROOT_INTEGRATION SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_ids=(ctest_it_1t ctest_it_2t)
# Analysis id reserved for this test.
readonly analysis_id=ctest_io_root

if [[ $# -ne 2 ]]; then
    echo "usage: $0 IO_ROOT_INTEGRATION SOURCE_DIR" >&2
    exit 2
fi
test_exe="$1"
source_dir="$2"

root="${RIPTIDE_DATA_DIR:-}"
[[ -n "${root}" ]] || root="${default_root}"

if ! "${source_dir}/scripts/data_root.sh" check --data-dir "${root}" --need 1; then
    echo "SKIP: data root ${root} is not usable, see the message above" >&2
    exit "${skip_code}"
fi
for id in "${run_ids[@]}"; do
    if [[ ! -f "${root}/sim/${id}/run.json" ]]; then
        echo "SKIP: run ${id} not found; it is written by the sim_integration test" >&2
        exit "${skip_code}"
    fi
done

scratch="${root}/analysis/${run_ids[0]}/${analysis_id}"
# Guard the removal: only an absolute path ending in the reserved ids.
if [[ "${scratch}" != /*"/analysis/${run_ids[0]}/${analysis_id}" ]]; then
    echo "refusing to remove unexpected path ${scratch}" >&2
    exit 1
fi
rm -rf -- "${scratch}"
mkdir -p -- "${scratch}"

RIPTIDE_IT_DATA_ROOT="${root}" RIPTIDE_IT_SCRATCH="${scratch}" \
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${test_exe}"
