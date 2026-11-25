#!/usr/bin/env bash
# Smoke test of riptide_case_b_geant4 on the 1-thread run of the riptide_sim
# integration test (ctest fixture sim_integration, a pencil beam along +Z):
# the same settings must give the same bytes, whatever the thread count; the
# toy source must run with the same protocol; --acceptance must be recorded
# and a missing acceptance.json refused; an existing stats directory must be
# refused and two sources must be a usage error. Its outputs go under
# reserved ids that are removed first.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the integration run.
#
# Usage: run_case_b_it.sh RIPTIDE_CASE_B_GEANT4 SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
readonly cuts=single_np,has_proton,contained
readonly direction=0,0,1
# Stats ids reserved for this test; the last one is never written.
readonly stats_ids=(ctest_case_b_a ctest_case_b_b ctest_case_b_toy ctest_case_b_acc_toy
    ctest_case_b_acc_g4 ctest_case_b_none)
readonly fixture_id=ctest_case_b_acc_fixture

if [[ $# -ne 2 ]]; then
    echo "usage: $0 RIPTIDE_CASE_B_GEANT4 SOURCE_DIR" >&2
    exit 2
fi
case_b_exe="$1"
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

remove_reserved() { # dir suffix
    # Guard the removal: only an absolute path ending in a reserved id.
    if [[ "$1" != /*"$2" ]]; then
        echo "refusing to remove unexpected path $1" >&2
        exit 1
    fi
    rm -rf -- "$1"
}
for id in "${stats_ids[@]}" "${fixture_id}"; do
    remove_reserved "${root}/stats/${id}" "/stats/${id}"
done

case_b() { # stats_id, then the source options
    local id="$1"
    shift
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${case_b_exe}" --data-dir "${root}" \
        --stats-id "${id}" --pilot "$@"
}
case_file() { echo "${root}/stats/$1/case_b.json"; }
# The records differ only in their stats id.
normalised() { # file id
    sed "s/\"$2\"/\"STATS_ID\"/" "$1"
}
expect_status() { # expected status, message, command...
    local expected="$1" message="$2"
    shift 2
    local status=0
    "$@" 2>/dev/null || status=$?
    if [[ "${status}" -ne "${expected}" ]]; then
        echo "FAIL: ${message} (exit ${status}, expected ${expected})" >&2
        exit 1
    fi
}

case_b "${stats_ids[0]}" --run-id "${run_id}" --threads 1
case_b "${stats_ids[1]}" --run-id "${run_id}" --threads 3
case_b "${stats_ids[2]}" --toy --direction "${direction}" --truth flat:10:19.5 --events 300

if ! cmp <(normalised "$(case_file "${stats_ids[0]}")" "${stats_ids[0]}") \
    <(normalised "$(case_file "${stats_ids[1]}")" "${stats_ids[1]}"); then
    echo "FAIL: the same settings on 1 and 3 threads gave different records" >&2
    exit 1
fi

expect_status 1 "an existing stats directory was not refused" \
    case_b "${stats_ids[0]}" --run-id "${run_id}"
expect_status 2 "two sources were not refused as a usage error" \
    case_b "${stats_ids[5]}" --run-id "${run_id}" --toy --direction "${direction}" \
    --truth flat:10:19.5 --events 10
expect_status 2 "a mono toy truth was not refused as a usage error" \
    case_b "${stats_ids[5]}" --toy --direction "${direction}" --truth mono:14 --events 10
expect_status 2 "a non-unit direction was not refused as a usage error" \
    case_b "${stats_ids[5]}" --toy --direction 0,0,2 --truth flat:10:19.5 --events 10
expect_status 1 "a missing acceptance.json was not an error" \
    case_b "${stats_ids[5]}" --toy --direction "${direction}" --truth flat:10:19.5 \
    --events 10 --acceptance "${stats_ids[5]}"
if [[ -e "${root}/stats/${stats_ids[5]}" ]]; then
    echo "FAIL: a refused call wrote its stats directory" >&2
    exit 1
fi

# A smooth acceptance table with no empty bin on the E_n domain of Case B.
fixture_dir="${root}/stats/${fixture_id}"
mkdir -p -- "${fixture_dir}"
python3 -I - "${fixture_dir}/acceptance.json" "${fixture_id}" <<'EOF'
import json
import math
import sys

path, stats_id = sys.argv[1:3]
En_edges = [10.0 + k for k in range(11)]
theta_edges = [(math.pi / 2) * t / 18 for t in range(18)] + [math.pi / 2]
counts = [float(1000 + 50 * t + 10 * k) for k in range(10) for t in range(18)]
record = {
    "acceptance_version": 1,
    "stats_id": stats_id,
    "source": {"run_id": "ctest_fixture"},
    "half": "all",
    "entries": sum(counts),
    "En_edges_MeV": En_edges,
    "theta_edges_rad": theta_edges,
    "counts": counts,
}
with open(path, "w") as f:
    json.dump(record, f, indent=2)
EOF
case_b "${stats_ids[3]}" --toy --direction "${direction}" --truth flat:10:19.5 --events 300 \
    --acceptance "${fixture_id}"
case_b "${stats_ids[4]}" --run-id "${run_id}" --acceptance "${fixture_id}"

python3 -I - "$(case_file "${stats_ids[0]}")" "$(case_file "${stats_ids[2]}")" \
    "$(case_file "${stats_ids[3]}")" "$(case_file "${stats_ids[4]}")" "${cuts}" \
    "${fixture_id}" <<'EOF'
import json
import sys

paths = sys.argv[1:5]
cuts, fixture_id = sys.argv[5:7]
geant4, toy, acc_toy, acc_geant4 = (json.load(open(p)) for p in paths)


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


top = {"case_b_version", "stats_id", "project_version", "git_describe", "source", "truth",
       "settings", "rows", "repeat", "contraction"}
records = (("geant4", geant4), ("toy", toy), ("acc_toy", acc_toy), ("acc_geant4", acc_geant4))
for name, record in records:
    if set(record) != top:
        fail(f"{name}: top-level keys {sorted(record)}")
    if record["case_b_version"] != 1:
        fail(f"{name}: version {record['case_b_version']}")
    if record["truth"]["omega_n"] != [0.0, 0.0, 1.0]:
        fail(f"{name}: omega_n {record['truth']['omega_n']}")
    if [row["events"] for row in record["rows"]] != [10, 30, 100, 300, 1000]:
        fail(f"{name}: rows {[row['events'] for row in record['rows']]}")
    if record["repeat"]["events"] != 1000:
        fail(f"{name}: repeat N = {record['repeat']['events']}")
    for row in record["rows"] + [record["repeat"]]:
        skipped = row["skipped"] is not None
        if skipped != ("ensemble" not in row):
            fail(f"{name}: row N = {row['events']} is half skipped")
        if skipped and row["experiments"] != 0:
            fail(f"{name}: skipped row with experiments")
        if not skipped:
            if not 2 <= row["experiments"] <= row["experiments_planned"]:
                fail(f"{name}: row N = {row['events']} has M = {row['experiments']}")
            if row["ensemble"]["entries"] != row["experiments"]:
                fail(f"{name}: row N = {row['events']} has entries != M")
            if row["ensemble"]["truth_in_cap_fraction"] != 1.0:
                fail(f"{name}: row N = {row['events']} lost the truth outside the cap")
    # The pool of the integration run is small: the repeat must be skipped,
    # and with it the slope that includes it.
    if record["repeat"]["skipped"] != "too few events":
        fail(f"{name}: the repeat was not skipped")
    if record["contraction"]["slope_with_repeat"] is not None:
        fail(f"{name}: a slope with a skipped repeat")

for name, record in (("geant4", geant4), ("toy", toy)):
    if record["settings"]["acceptance"] is not None:
        fail(f"{name}: acceptance {record['settings']['acceptance']}")
for name, record in (("acc_toy", acc_toy), ("acc_geant4", acc_geant4)):
    if record["settings"]["acceptance"]["stats_id"] != fixture_id:
        fail(f"{name}: acceptance {record['settings']['acceptance']}")

for name, record in (("geant4", geant4), ("acc_geant4", acc_geant4)):
    source = record["source"]
    if source["kind"] != "geant4" or source["cuts"] != cuts:
        fail(f"{name}: source {source['kind']}, cuts {source['cuts']}")
    if source["events_selected"] != source["cut_flow"][-1]["passed"]:
        fail(f"{name}: events_selected differs from the last cut")
    if record["truth"]["spectrum"] != {"mode": "mono", "En_MeV": 14.0}:
        fail(f"{name}: spectrum {record['truth']['spectrum']}")
if geant4["source"] != acc_geant4["source"]:
    fail("the acceptance changed the Geant4 pool")

for name, record in (("toy", toy), ("acc_toy", acc_toy)):
    if record["source"] != {"kind": "toy", "events": 300}:
        fail(f"{name}: source {record['source']}")
    if record["truth"]["spectrum"] != {"mode": "flat", "En_min_MeV": 10.0, "En_max_MeV": 19.5}:
        fail(f"{name}: spectrum {record['truth']['spectrum']}")
EOF
echo "riptide_case_b_geant4 smoke test passed"
