#!/usr/bin/env bash
# Smoke test of riptide_acceptance and of riptide_case_a_geant4 --acceptance
# (D40) on the 1-thread run of the riptide_sim integration test (ctest fixture
# sim_integration).
#
# riptide_acceptance: the same run and settings must give the same bytes, the
# entries must be the selected events of the Case A chain, the even and odd
# halves must add up to the whole, and an existing stats directory and a bad
# --half must be refused. The integration run is small (about a hundred
# selected events), so its table has empty bins: the Case A part uses instead
# a smooth table written here, with no empty bin, on the 10-20 MeV range.
#
# Its outputs go under reserved ids that are removed first. Skipped (exit 4,
# ctest SKIP_RETURN_CODE) without a usable data root or without the
# integration run.
#
# Usage: run_acceptance_it.sh RIPTIDE_ACCEPTANCE RIPTIDE_CASE_A_GEANT4 SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
readonly cuts=single_np,has_proton,contained
# Stats ids reserved for this test.
readonly acc_ids=(ctest_acc_a ctest_acc_b ctest_acc_even ctest_acc_odd ctest_acc_bad)
readonly fixture_id=ctest_acc_fixture
readonly case_ids=(ctest_acc_case_toy_a ctest_acc_case_toy_b ctest_acc_case_g4 ctest_acc_case_missing)

if [[ $# -ne 3 ]]; then
    echo "usage: $0 RIPTIDE_ACCEPTANCE RIPTIDE_CASE_A_GEANT4 SOURCE_DIR" >&2
    exit 2
fi
acceptance_exe="$1"
case_a_exe="$2"
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
for id in "${acc_ids[@]}" "${fixture_id}" "${case_ids[@]}"; do
    remove_reserved "${root}/stats/${id}" "/stats/${id}"
done

capped() {
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "$@"
}
acceptance() { # stats_id, then options
    local id="$1"
    shift
    capped "${acceptance_exe}" --data-dir "${root}" --stats-id "${id}" --run-id "${run_id}" "$@"
}
case_a() { # stats_id, then options
    local id="$1"
    shift
    capped "${case_a_exe}" --data-dir "${root}" --stats-id "${id}" --pilot "$@"
}
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

# --- riptide_acceptance ---------------------------------------------------
acceptance "${acc_ids[0]}"
acceptance "${acc_ids[1]}"
acceptance "${acc_ids[2]}" --half even
acceptance "${acc_ids[3]}" --half odd
acc_file() { echo "${root}/stats/$1/acceptance.json"; }

if ! cmp <(normalised "$(acc_file "${acc_ids[0]}")" "${acc_ids[0]}") \
    <(normalised "$(acc_file "${acc_ids[1]}")" "${acc_ids[1]}"); then
    echo "FAIL: the same run and settings gave different acceptance records" >&2
    exit 1
fi
expect_status 1 "an existing stats directory was not refused" acceptance "${acc_ids[0]}"
expect_status 2 "a bad --half was not a usage error" acceptance "${acc_ids[4]}" --half half
if [[ -e "${root}/stats/${acc_ids[4]}" ]]; then
    echo "FAIL: a usage error wrote ${acc_ids[4]}" >&2
    exit 1
fi

python3 -I - "$(acc_file "${acc_ids[0]}")" "$(acc_file "${acc_ids[2]}")" \
    "$(acc_file "${acc_ids[3]}")" "${run_id}" "${cuts}" <<'EOF'
import json
import math
import sys

all_path, even_path, odd_path, run_id, cuts = sys.argv[1:6]
records = {}
for name, path in (("all", all_path), ("even", even_path), ("odd", odd_path)):
    with open(path) as f:
        records[name] = json.load(f)


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


top = {"acceptance_version", "stats_id", "project_version", "git_describe", "source", "cuts",
       "cut_flow", "events_read", "events_selected", "half", "entries", "En_edges_MeV",
       "theta_edges_rad", "counts"}
for name, record in records.items():
    if set(record) != top:
        fail(f"{name}: top-level keys {sorted(record)}")
    if record["half"] != name or record["cuts"] != cuts or record["source"]["run_id"] != run_id:
        fail(f"{name}: half {record['half']}, cuts {record['cuts']}")
    if record["events_selected"] != record["cut_flow"][-1]["passed"]:
        fail(f"{name}: events_selected differs from the last cut")
    if record["entries"] != sum(record["counts"]):
        fail(f"{name}: entries differ from the sum of the counts")
    theta = record["theta_edges_rad"]
    if len(theta) != 19 or theta[0] != 0.0 or not math.isclose(theta[-1], math.pi / 2, rel_tol=1e-15):
        fail(f"{name}: theta_p edges {theta}")
    # The integration run is mono at 14 MeV: one bin centred on it.
    if record["En_edges_MeV"] != [13.5, 14.5]:
        fail(f"{name}: E_n edges {record['En_edges_MeV']}")
    if len(record["counts"]) != 18:
        fail(f"{name}: {len(record['counts'])} counts")

whole = records["all"]
if whole["entries"] != whole["events_selected"] or whole["entries"] <= 0:
    fail(f"all: entries {whole['entries']}, selected {whole['events_selected']}")
for i, count in enumerate(whole["counts"]):
    if records["even"]["counts"][i] + records["odd"]["counts"][i] != count:
        fail(f"the halves do not add up in bin {i}")
EOF

# --- riptide_case_a_geant4 --acceptance -------------------------------------
# A smooth table with no empty bin: the forward deficit of D40 in shape, on
# the E_n grid range of Case A.
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

case_a "${case_ids[0]}" --toy --truth mono:14 --events 300 --acceptance "${fixture_id}"
case_a "${case_ids[1]}" --toy --truth mono:14 --events 300 --acceptance "${fixture_id}"
case_a "${case_ids[2]}" --run-id "${run_id}" --acceptance "${fixture_id}"
case_file() { echo "${root}/stats/$1/case_a.json"; }
if ! cmp <(normalised "$(case_file "${case_ids[0]}")" "${case_ids[0]}") \
    <(normalised "$(case_file "${case_ids[1]}")" "${case_ids[1]}"); then
    echo "FAIL: the same settings with --acceptance gave different records" >&2
    exit 1
fi
expect_status 1 "a missing acceptance.json was not an error" \
    case_a "${case_ids[3]}" --toy --truth mono:14 --events 10 --acceptance "${acc_ids[4]}"

python3 -I - "$(case_file "${case_ids[0]}")" "$(case_file "${case_ids[2]}")" "${fixture_id}" <<'EOF'
import json
import sys

toy_path, geant4_path, fixture_id = sys.argv[1:4]


def fail(message):
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


for name, path in (("toy", toy_path), ("geant4", geant4_path)):
    with open(path) as f:
        record = json.load(f)
    if record["case_a_version"] != 2:
        fail(f"{name}: case_a_version {record['case_a_version']}")
    acc = record["settings"]["acceptance"]
    expected = {"stats_id", "acceptance_version", "run_id", "half", "entries", "En_bins",
                "theta_bins", "ratio_max"}
    if acc is None or set(acc) != expected:
        fail(f"{name}: settings.acceptance {acc}")
    if acc["stats_id"] != fixture_id or acc["run_id"] != "ctest_fixture":
        fail(f"{name}: acceptance from {acc['stats_id']}, run {acc['run_id']}")
    if acc["En_bins"] != 10 or acc["theta_bins"] != 18 or not acc["ratio_max"] > 1.0:
        fail(f"{name}: acceptance bins {acc['En_bins']} x {acc['theta_bins']}, max {acc['ratio_max']}")
EOF
echo "riptide_acceptance smoke test passed"
