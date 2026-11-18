#!/usr/bin/env bash
# Smoke test of riptide_sim_check on the 1-thread run of the integration test
# (ctest fixture sim_integration): the app must finish, write every output
# file and report zero consistency violations. Physics numbers are not judged
# here; this checks that the pipeline holds together.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the evaluated reference files.
#
# Usage: run_check_smoke.sh RIPTIDE_SIM_CHECK SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
# Check id reserved for this test; its directory is recreated every time.
readonly check_id=ctest_smoke
readonly reference_files=(
    jeff33/n_0125_1-H-1.dat jeff33/n_0600_6-C-0.dat jeff33/n_0628_6-C-13.dat
    endfb80/n-001_H_001.endf endfb80/n-006_C_012.endf endfb80/n-006_C_013.endf)

if [[ $# -ne 2 ]]; then
    echo "usage: $0 RIPTIDE_SIM_CHECK SOURCE_DIR" >&2
    exit 2
fi
check_exe="$1"
source_dir="$2"

root="${RIPTIDE_DATA_DIR:-}"
[[ -n "${root}" ]] || root="${default_root}"

if ! "${source_dir}/scripts/data_root.sh" check --data-dir "${root}" --need 1; then
    echo "SKIP: data root ${root} is not usable, see the message above" >&2
    exit "${skip_code}"
fi
for f in "${reference_files[@]}"; do
    if [[ ! -f "${root}/reference/${f}" ]]; then
        echo "SKIP: reference file ${root}/reference/${f} not found" >&2
        exit "${skip_code}"
    fi
done

out_dir="${root}/validation/${run_id}/${check_id}"
# Guard the removal: only an absolute path ending in the reserved ids.
if [[ "${out_dir}" != /*"/validation/${run_id}/${check_id}" ]]; then
    echo "refusing to remove unexpected path ${out_dir}" >&2
    exit 1
fi
rm -rf -- "${out_dir}"

"${source_dir}/scripts/run_capped.sh" -m 5G -- "${check_exe}" \
    --data-dir "${root}" --run-id "${run_id}" --check-id "${check_id}"

for f in checks.json spectrum_u.csv depth_z.csv fractions.csv legendre.csv residuals_Ep.csv; do
    [[ -s "${out_dir}/${f}" ]] || { echo "missing or empty ${out_dir}/${f}" >&2; exit 1; }
done

# The sha256 recorded for the first reference file must match sha256sum's.
python3 -I - "${out_dir}/checks.json" "${root}/reference/${reference_files[0]}" \
    "$(sha256sum "${root}/reference/${reference_files[0]}" | cut -d' ' -f1)" <<'PY'
import json
import sys

path, ref_path, ref_sha = sys.argv[1:4]
with open(path, encoding="utf-8") as f:
    doc = json.load(f)
for key in ("checks_version", "run_id", "check_id", "versions", "reference_files",
            "hp_elastic_H1", "medium", "beam", "results"):
    assert key in doc, f"checks.json has no key {key}"
assert doc["checks_version"] == 3, f"unexpected checks_version {doc['checks_version']}"
results = doc["results"]
for key in ("spectrum_u", "depth_z", "fractions", "legendre", "Ep_residuals",
            "Ep_physical_offset", "counters"):
    assert key in results, f"results has no key {key}"
# The HP target mass of H-1 as G4NDL 4.7.1 gives it, in neutron masses.
hp = doc["hp_elastic_H1"]
assert abs(hp["awr"] - 0.999167) < 1e-12, f"unexpected H-1 AWR {hp['awr']}"
# Oracle: with that mass the simulated Ep is reproduced to rounding. The
# bound is relative to the largest proton energy of the source, a slight
# overestimate of Ep_max(En_max) that does not matter at this tolerance.
En_max = doc["beam"]["En_max_MeV"]
residual = results["Ep_residuals"]["max_abs"]
assert residual < 1e-9 * En_max, f"Ep residual {residual} MeV is not rounding"
# The departure from the bare-proton kinematics is the HP convention, ~m_e/m_p.
offset = results["Ep_physical_offset"]["relative_mean"]
assert 5.3e-4 <= offset <= 5.5e-4, f"relative Ep offset {offset} outside [5.3e-4, 5.5e-4]"
total = results["counters"]["must_be_zero"]["total"]
assert total == 0, f"{total} consistency violations: {results['counters']['must_be_zero']}"
assert results["events"] > 0, "no events read"
# JEFF-3.3 is composed as the HP model composes carbon: C-0 for C-12, plus C-13.
jeff = [r["nuclide"] for r in doc["reference_files"] if r["library"] == "jeff33"]
assert jeff == ["H-1", "C-0", "C-13"], f"unexpected JEFF-3.3 files {jeff}"
recorded = {r["path"]: r["sha256"] for r in doc["reference_files"]}
assert recorded.get(ref_path) == ref_sha, f"sha256 of {ref_path}: {recorded.get(ref_path)} != {ref_sha}"
print(f"checks.json ok: {results['events']} events, {results['np_events']} n-p")
PY
