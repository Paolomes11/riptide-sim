#!/usr/bin/env bash
# Integration test of riptide_sim: three short runs of the same events, with 1
# and 2 worker threads and, for the kB = 0 oracle, with 1 thread and the Birks
# constant set to 0; then a Catch2 reader that checks the written files.
#
# The runs go to the data root (the external SSD): the system disk does not
# have the free-space margin that riptide_sim requires. Without a usable data
# root the test is skipped (exit 4, ctest SKIP_RETURN_CODE), not failed.
#
# Usage: run_integration.sh RIPTIDE_SIM READER SOURCE_DIR
#   The data root is resolved like the simulation does: $RIPTIDE_DATA_DIR,
#   then /mnt/external_ssd/riptide-sim.
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly events=1000
readonly seed=20260907
# Run ids reserved for this test; their directories are recreated every time.
readonly run_ids=(ctest_it_1t ctest_it_2t ctest_it_kB0)
readonly threads=(1 2 1)

if [[ $# -ne 3 ]]; then
    echo "usage: $0 RIPTIDE_SIM READER SOURCE_DIR" >&2
    exit 2
fi
sim_exe="$1"
reader_exe="$2"
source_dir="$3"

root="${RIPTIDE_DATA_DIR:-}"
[[ -n "${root}" ]] || root="${default_root}"

if ! "${source_dir}/scripts/data_root.sh" check --data-dir "${root}" --need 1; then
    echo "SKIP: data root ${root} is not usable, see the message above" >&2
    exit "${skip_code}"
fi

# The kB = 0 config is mono_14MeV.json with only the Birks constant changed,
# written outside the repository and removed on exit.
config_dir="$(mktemp -d)"
trap 'rm -rf -- "${config_dir}"' EXIT
kB0_config="${config_dir}/mono_14MeV_kB0.json"
python3 -I - "${source_dir}/config/mono_14MeV.json" "${kB0_config}" <<'PY'
import json
import sys

src, dst = sys.argv[1:3]
with open(src, encoding="utf-8") as f:
    config = json.load(f)
config["scintillator"] = {"birks_kB_mm_per_MeV": 0.0,
                          "birks_source": "integration test: kB = 0 oracle"}
with open(dst, "w", encoding="utf-8") as f:
    json.dump(config, f, indent=2)
PY
configs=("${source_dir}/config/mono_14MeV.json" "${source_dir}/config/mono_14MeV.json"
    "${kB0_config}")

run_dirs=()
for i in "${!run_ids[@]}"; do
    run_dir="${root}/sim/${run_ids[$i]}"
    # Guard the removal: only an absolute path ending in the reserved run id.
    if [[ "${run_dir}" != /*"/sim/${run_ids[$i]}" ]]; then
        echo "refusing to remove unexpected path ${run_dir}" >&2
        exit 1
    fi
    rm -rf -- "${run_dir}"
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${sim_exe}" \
        --config "${configs[$i]}" \
        --data-dir "${root}" --run-id "${run_ids[$i]}" \
        --events "${events}" --seed "${seed}" --threads "${threads[$i]}"
    run_dirs+=("${run_dir}")
done

RIPTIDE_IT_RUN_1T="${run_dirs[0]}" RIPTIDE_IT_RUN_2T="${run_dirs[1]}" \
    RIPTIDE_IT_RUN_KB0="${run_dirs[2]}" \
    RIPTIDE_IT_EVENTS="${events}" "${reader_exe}"
