#!/usr/bin/env bash
# Smoke test of riptide_case_a_geant4 on the 1-thread run of the riptide_sim
# integration test (ctest fixture sim_integration): the same settings must
# give the same bytes, the toy source must run with the same protocol, an
# existing stats directory must be refused, and the cut flow must be the one
# riptide_analysis_prepare records for the same chain. Its outputs go under
# reserved ids that are removed first.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the integration run.
#
# Usage: run_case_a_it.sh RIPTIDE_CASE_A_GEANT4 RIPTIDE_ANALYSIS_PREPARE SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
readonly cuts=single_np,has_proton,contained
# Stats ids and the analysis id reserved for this test.
readonly stats_ids=(ctest_case_a_a ctest_case_a_b ctest_case_a_c ctest_case_a_d)
readonly analysis_id=ctest_case_a_prepare

if [[ $# -ne 3 ]]; then
    echo "usage: $0 RIPTIDE_CASE_A_GEANT4 RIPTIDE_ANALYSIS_PREPARE SOURCE_DIR" >&2
    exit 2
fi
case_a_exe="$1"
prepare_exe="$2"
source_dir="$3"

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

remove_reserved() { # dir suffix
    # Guard the removal: only an absolute path ending in a reserved id.
    if [[ "$1" != /*"$2" ]]; then
        echo "refusing to remove unexpected path $1" >&2
        exit 1
    fi
    rm -rf -- "$1"
}
for id in "${stats_ids[@]}"; do
    remove_reserved "${root}/stats/${id}" "/stats/${id}"
done
remove_reserved "${root}/analysis/${run_id}/${analysis_id}" "/analysis/${run_id}/${analysis_id}"

case_a() { # stats_id, then the source options
    local id="$1"
    shift
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${case_a_exe}" --data-dir "${root}" \
        --stats-id "${id}" --pilot "$@"
}

case_a "${stats_ids[0]}" --run-id "${run_id}"
case_a "${stats_ids[1]}" --run-id "${run_id}"
case_a "${stats_ids[2]}" --toy --truth mono:14 --events 300

file_a="${root}/stats/${stats_ids[0]}/case_a.json"
file_b="${root}/stats/${stats_ids[1]}/case_a.json"
file_toy="${root}/stats/${stats_ids[2]}/case_a.json"

# The records differ only in their stats id.
normalised() { # file id
    sed "s/\"$2\"/\"STATS_ID\"/" "$1"
}
if ! cmp <(normalised "${file_a}" "${stats_ids[0]}") <(normalised "${file_b}" "${stats_ids[1]}"); then
    echo "FAIL: the same settings gave different records" >&2
    exit 1
fi

# An existing stats directory is refused (runtime error, exit 1).
status=0
case_a "${stats_ids[0]}" --run-id "${run_id}" 2>/dev/null || status=$?
if [[ "${status}" -ne 1 ]]; then
    echo "FAIL: an existing stats directory was not refused (exit ${status})" >&2
    exit 1
fi
# Two sources at once are a usage error (exit 2), and nothing is written.
status=0
case_a "${stats_ids[3]}" --run-id "${run_id}" --toy --truth mono:14 --events 10 \
    2>/dev/null || status=$?
if [[ "${status}" -ne 2 || -e "${root}/stats/${stats_ids[3]}" ]]; then
    echo "FAIL: two sources were not refused as a usage error (exit ${status})" >&2
    exit 1
fi

# The same chain in riptide_analysis_prepare gives the cut flow to compare.
"${source_dir}/scripts/run_capped.sh" -m 5G -- "${prepare_exe}" --data-dir "${root}" \
    --run-id "${run_id}" --analysis-id "${analysis_id}" --cuts "${cuts}" \
    --seed 20260907 --sigma-Ep 0.10 --sigma-theta 0.08
prepare_record="${root}/analysis/${run_id}/${analysis_id}/analysis.json"

python3 -I - "${file_a}" "${file_toy}" "${prepare_record}" "${cuts}" <<'EOF'
import json
import sys

geant4_path, toy_path, prepare_path, cuts = sys.argv[1:5]
with open(geant4_path) as f:
    geant4 = json.load(f)
with open(toy_path) as f:
    toy = json.load(f)
with open(prepare_path) as f:
    prepare = json.load(f)["selection"]


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


top = {"case_a_version", "stats_id", "project_version", "git_describe", "source", "truth",
       "settings", "contraction", "theta_prior", "energy_prior"}
for name, record in (("geant4", geant4), ("toy", toy)):
    if set(record) != top:
        fail(f"{name}: top-level keys {sorted(record)}")
    if record["case_a_version"] != 2 or record["settings"]["acceptance"] is not None:
        fail(f"{name}: version {record['case_a_version']}, acceptance {record['settings']['acceptance']}")
    if record["truth"] != {"mode": "mono", "En_MeV": 14.0}:
        fail(f"{name}: truth {record['truth']}")
    for section in ("contraction", "energy_prior"):
        for row in record[section]["rows"]:
            skipped = row["skipped"] is not None
            if skipped != ("ensemble" not in row and "reference" not in row):
                fail(f"{name}: {section} row N = {row['events']} is half skipped")
            if skipped and row["experiments"] != 0:
                fail(f"{name}: skipped {section} row with experiments")
            if not skipped and not 2 <= row["experiments"] <= row["experiments_planned"]:
                fail(f"{name}: {section} row N = {row['events']} has M = {row['experiments']}")

source = geant4["source"]
if source["kind"] != "geant4" or source["cuts"] != cuts:
    fail(f"geant4 source {source['kind']}, cuts {source['cuts']}")
# The pool of the integration run is small: the largest rows must be skipped,
# and that must be visible in the file.
if geant4["theta_prior"]["skipped"] != "too few events":
    fail("geant4: the theta_p prior row with N = 1000 was not skipped")
if source["events_selected"] != source["cut_flow"][-1]["passed"]:
    fail("geant4: events_selected differs from the last cut")
if source["cut_flow"] != prepare["cuts"] or source["cuts"] != prepare["text"] or source["events_read"] != prepare["events_read"]:
    fail("geant4: cut flow differs from riptide_analysis_prepare")
if source["events_selected"] != prepare["events_selected"]:
    fail("geant4: selected events differ from riptide_analysis_prepare")

if toy["source"] != {"kind": "toy", "events": 300}:
    fail(f"toy source {toy['source']}")
EOF
echo "riptide_case_a_geant4 smoke test passed"
