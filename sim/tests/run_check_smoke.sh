#!/usr/bin/env bash
# Smoke test of riptide_sim_check on the 1-thread run of the integration test
# (ctest fixture sim_integration): the app must finish, write every output
# file and report zero consistency violations. Physics numbers are not judged
# here; this checks that the pipeline holds together.
#
# The checker then reads the 2-thread run of the same events: its outputs must
# be byte-identical to the 1-thread ones, ids aside. The rows are merged by
# event_id, so neither the number of files nor how Geant4 spread the events
# over them may change a sum.
#
# Skipped (exit 4, ctest SKIP_RETURN_CODE) without a usable data root or
# without the evaluated reference files.
#
# Usage: run_check_smoke.sh RIPTIDE_SIM_CHECK SOURCE_DIR
set -euo pipefail

readonly skip_code=4
readonly default_root="/mnt/external_ssd/riptide-sim"
readonly run_id=ctest_it_1t
readonly run_id_2t=ctest_it_2t
# Check id reserved for this test; its directories are recreated every time.
readonly check_id=ctest_smoke
readonly output_files=(checks.json spectrum_u.csv depth_z.csv fractions.csv legendre.csv
    model_mix.csv residuals_Ep.csv deposit_spectra.csv proton_exits.csv pstar_comparison.csv)
readonly reference_files=(
    jeff33/n_0125_1-H-1.dat jeff33/n_0600_6-C-0.dat jeff33/n_0628_6-C-13.dat
    endfb80/n-001_H_001.endf endfb80/n-006_C_012.endf endfb80/n-006_C_013.endf
    pstar/pstar_216.txt)
readonly pstar_file="pstar/pstar_216.txt"

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

# Runs the checker on run $1 into a fresh validation/$1/${check_id}.
run_check() {
    local dir="${root}/validation/$1/${check_id}"
    # Guard the removal: only an absolute path ending in the reserved ids.
    if [[ "${dir}" != /*"/validation/$1/${check_id}" ]]; then
        echo "refusing to remove unexpected path ${dir}" >&2
        exit 1
    fi
    rm -rf -- "${dir}"
    "${source_dir}/scripts/run_capped.sh" -m 5G -- "${check_exe}" \
        --data-dir "${root}" --run-id "$1" --check-id "${check_id}"
    for f in "${output_files[@]}"; do
        [[ -s "${dir}/${f}" ]] || { echo "missing or empty ${dir}/${f}" >&2; exit 1; }
    done
}

out_dir="${root}/validation/${run_id}/${check_id}"
out_dir_2t="${root}/validation/${run_id_2t}/${check_id}"
run_check "${run_id}"
run_check "${run_id_2t}"

for f in "${output_files[@]}"; do
    [[ "${f}" == checks.json ]] && continue
    cmp -- "${out_dir}/${f}" "${out_dir_2t}/${f}" ||
        { echo "${f} differs between ${run_id} and ${run_id_2t}" >&2; exit 1; }
done
# checks.json names its run; everything else must match exactly.
python3 -I - "${out_dir}/checks.json" "${out_dir_2t}/checks.json" <<'PY'
import json
import sys

docs = []
for path in sys.argv[1:3]:
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    doc.pop("run_id")
    docs.append(doc)
assert docs[0] == docs[1], "checks.json differs between the 1- and 2-thread runs"
PY

# The sha256 recorded for the first reference file and the PSTAR table must
# match sha256sum's.
python3 -I - "${out_dir}/checks.json" "${root}/reference/${reference_files[0]}" \
    "$(sha256sum "${root}/reference/${reference_files[0]}" | cut -d' ' -f1)" \
    "${root}/reference/${pstar_file}" \
    "$(sha256sum "${root}/reference/${pstar_file}" | cut -d' ' -f1)" <<'PY'
import json
import sys

path, ref_path, ref_sha, pstar_path, pstar_sha = sys.argv[1:6]
with open(path, encoding="utf-8") as f:
    doc = json.load(f)
for key in ("checks_version", "run_id", "check_id", "versions", "reference_files",
            "hp_elastic_H1", "hp_deuterium", "medium", "beam", "results"):
    assert key in doc, f"checks.json has no key {key}"
assert doc["checks_version"] == 9, f"unexpected checks_version {doc['checks_version']}"
results = doc["results"]
for key in ("spectrum_u", "depth_z", "fractions", "legendre", "Ep_residuals",
            "Ep_physical_offset", "model_mix", "counters", "deposits"):
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
# A 14 MeV beam lies below the HP/CHIPS window: every n-p event is HP.
mix = results["model_mix"]
assert mix["observed_non_hp"] == 0, f"{mix['observed_non_hp']} non-HP n-p events at 14 MeV"
assert results["np_hp_events"] == results["np_events"], "n-p events left out of the HP sample"
total = results["counters"]["must_be_zero"]["total"]
assert total == 0, f"{total} consistency violations: {results['counters']['must_be_zero']}"
assert results["events"] > 0, "no events read"
# JEFF-3.3 is composed as the HP model composes carbon: C-0 for C-12, plus C-13.
jeff = [r["nuclide"] for r in doc["reference_files"] if r["library"] == "jeff33"]
assert jeff == ["H-1", "C-0", "C-13"], f"unexpected JEFF-3.3 files {jeff}"
recorded = {r["path"]: r["sha256"] for r in doc["reference_files"]}
assert recorded.get(ref_path) == ref_sha, f"sha256 of {ref_path}: {recorded.get(ref_path)} != {ref_sha}"
# H-2 enters the expectations from its three G4NDL cross-section files.
h2 = [f["file"].split("/")[-3] for f in doc["hp_deuterium"]
      if f["file"].endswith("/CrossSection/1_2_Hydrogen.z") and len(f["sha256"]) == 64]
assert h2 == ["Elastic", "Inelastic", "Capture"], f"unexpected H-2 files {doc['hp_deuterium']}"
# The integration run is schema v3 with the nominal kB: the deposit checks run.
deposits = results["deposits"]
assert deposits["skipped"] is None, f"deposit checks skipped: {deposits['skipped']}"
assert deposits["events"] == results["events"], "deposit checks saw a different event count"
assert sum(deposits["event_classes"].values()) == deposits["events"], "event classes do not add up"
assert deposits["setup"]["birks_kB_mm_per_MeV"] > 0, "the integration run should use the nominal kB"
dep_total = deposits["counters"]["must_be_zero"]["total"]
assert dep_total == 0, f"{dep_total} deposit violations: {deposits['counters']['must_be_zero']}"
pstar = deposits["pstar"]
assert pstar["skipped"] is None, f"PSTAR comparison skipped: {pstar['skipped']}"
assert pstar["file"] == pstar_path, f"PSTAR table {pstar['file']} instead of {pstar_path}"
assert pstar["sha256"] == pstar_sha, "sha256 of the PSTAR table differs from sha256sum"
# The compared protons and those below the table are the contained balance sample.
contained = deposits["contained_balance_MeV"]["samples"]
assert pstar["samples"] + pstar["below_table"] == contained, \
    f"PSTAR saw {pstar['samples']} + {pstar['below_table']} protons, the balance {contained}"
assert pstar["range_mm"]["samples"] == pstar["samples"], "range residuals miss protons"
assert len(pstar["bins"]) == deposits["setup"]["bins"], "PSTAR bins differ from the deposit bins"
print(f"checks.json ok: {results['events']} events, {results['np_events']} n-p")
PY
