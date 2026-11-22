#!/usr/bin/env bash
# Smoke test of riptide_analysis_prepare on the 1-thread run of the riptide_sim
# integration test (ctest fixture sim_integration): the same seed must give
# the same record, a different seed a different one, and an existing analysis
# directory must be refused. Its outputs go under reserved analysis ids that
# are removed first.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the integration run.
#
# Usage: run_prepare_it.sh RIPTIDE_ANALYSIS_PREPARE SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
# Analysis ids reserved for this test.
readonly analysis_ids=(ctest_prepare_a ctest_prepare_b ctest_prepare_c)
readonly seed=20260907
readonly other_seed=20260908

if [[ $# -ne 2 ]]; then
    echo "usage: $0 RIPTIDE_ANALYSIS_PREPARE SOURCE_DIR" >&2
    exit 2
fi
prepare_exe="$1"
source_dir="$2"

root="${RIPTIDE_DATA_DIR:-}"
[[ -n "${root}" ]] || root="${default_root}"

if ! "${source_dir}/scripts/data_root.sh" check --data-dir "${root}" --need 1; then
    echo "SKIP: data root ${root} is not usable, see the message above" >&2
    exit "${skip_code}"
fi
if [[ ! -f "${root}/sim/${run_id}/run.json" ]]; then
    echo "SKIP: run ${run_id} not found; it is written by the sim_integration test" >&2
    exit "${skip_code}"
fi

for id in "${analysis_ids[@]}"; do
    dir="${root}/analysis/${run_id}/${id}"
    # Guard the removal: only an absolute path ending in the reserved ids.
    if [[ "${dir}" != /*"/analysis/${run_id}/${id}" ]]; then
        echo "refusing to remove unexpected path ${dir}" >&2
        exit 1
    fi
    rm -rf -- "${dir}"
done

prepare() { # analysis_id seed
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${prepare_exe}" --data-dir "${root}" \
        --run-id "${run_id}" --analysis-id "$1" --seed "$2" --sigma-Ep 0.10 --sigma-theta 0.08 \
        --cuts single_np,contained
}

prepare "${analysis_ids[0]}" "${seed}"
prepare "${analysis_ids[1]}" "${seed}"
prepare "${analysis_ids[2]}" "${other_seed}"

dir_a="${root}/analysis/${run_id}/${analysis_ids[0]}"
dir_b="${root}/analysis/${run_id}/${analysis_ids[1]}"
dir_c="${root}/analysis/${run_id}/${analysis_ids[2]}"
for key in '"hp_conversion"' '"resolution"' '"pull_variance"'; do
    if ! grep -q -- "${key}" "${dir_a}/analysis.json"; then
        echo "FAIL: ${key} missing from analysis.json" >&2
        exit 1
    fi
done

# The records differ only in their analysis id.
normalised() { # dir id
    sed "s/\"$2\"/\"ANALYSIS_ID\"/" "$1/analysis.json"
}
if ! diff <(normalised "${dir_a}" "${analysis_ids[0]}") \
        <(normalised "${dir_b}" "${analysis_ids[1]}"); then
    echo "FAIL: the same seed gave different records" >&2
    exit 1
fi
cmp "${dir_a}/cutflow.csv" "${dir_b}/cutflow.csv"
if diff -q <(normalised "${dir_a}" "${analysis_ids[0]}") \
        <(normalised "${dir_c}" "${analysis_ids[2]}") >/dev/null; then
    echo "FAIL: a different seed gave the same record" >&2
    exit 1
fi
# The cut flow does not depend on the seed.
cmp "${dir_a}/cutflow.csv" "${dir_c}/cutflow.csv"

# An existing analysis directory is refused (runtime error, exit 1).
status=0
prepare "${analysis_ids[0]}" "${seed}" 2>/dev/null || status=$?
if [[ "${status}" -ne 1 ]]; then
    echo "FAIL: an existing analysis directory was not refused (exit ${status})" >&2
    exit 1
fi
echo "riptide_analysis_prepare smoke test passed"
