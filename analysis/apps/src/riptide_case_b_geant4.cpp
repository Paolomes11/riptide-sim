// Case B on Geant4 events (D45): reads one riptide_sim run with a fixed
// neutron direction, applies the cut chain of Case A, converts the HP proton
// energies, draws the detector resolution with a recorded seed and splits the
// selected events into disjoint experiments, then runs the rows of the 7a
// checklist on them with the uniform direction prior. With --toy the pool
// comes from the Case B toy instead, with the same direction, truth and
// protocol, so the difference of the two outputs is the cost of the detector
// realism. With --acceptance the model carries the selection acceptance of a
// flat run (D40) and the toy draws its events with it, so it stays the
// perfect-model twin of Geant4. Experiments run on several threads (D44);
// the numbers do not depend on --threads. It writes case_b.json under the
// stats directory of the data root; the file holds no wall time, so two runs
// with the same settings give the same bytes. Timings go to stdout only.

#include <lyra/lyra.hpp>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "riptide/apps/case_a_geant4.hpp"
#include "riptide/apps/case_b_geant4.hpp"
#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/direction.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/random.hpp"
#include "riptide/common/units.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/acceptance_record.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"
#include "riptide/stats/acceptance.hpp"
#include "riptide/stats/calibration_b.hpp"
#include "riptide/stats/case_b.hpp"
#include "riptide/stats/case_b_posterior.hpp"
#include "riptide/stats/direction_grid.hpp"
#include "riptide/stats/grid.hpp"
#include "riptide/stats/priors.hpp"
#include "riptide/stats/pseudo_experiment_b.hpp"
#include "riptide/stats/split.hpp"

namespace {

namespace fs = std::filesystem;
namespace apps = riptide::apps;
namespace io = riptide::io;
namespace st = riptide::stats;
using riptide::geometry::Direction;
using Clock = std::chrono::steady_clock;
using nlohmann::json;

constexpr const char* kAppName = "riptide_case_b_geant4";
constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
constexpr std::uint64_t kProjectSeed = 20260907;
constexpr st::CaseBResolution kDefaultResolution{0.10, 0.08}; // MeV, rad
constexpr unsigned kDefaultThreads = 4;
constexpr int kCaseBVersion = 1;
constexpr const char* kCaseBName = "case_b.json";
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

// The cut chain of Case A (D30, D34): one n-p scattering, a proton, and the
// proton contained, so its deposit is its energy and its track is measured.
constexpr const char* kCaseBCuts = "single_np,has_proton,contained";

// The rows of the 7a checklist (7_predictions.md, section 4.1).
struct SizedRow {
    std::size_t events;      // N per experiment
    std::size_t experiments; // planned M
};
const std::vector<SizedRow> kChecklistRows{
    {10, 400}, {30, 300}, {100, 200}, {300, 100}, {1000, 40}};
// The repeat at the N of the last row, on indices disjoint from it.
constexpr std::size_t kRepeatExperiments = 240;
// The contraction slope is fitted from this N up, as in 7a.
constexpr std::size_t kSlopeMinEvents = 30;

// Sections of the split engine, one per row; the last row and the repeat
// share one permutation, so they share no event (D37).
constexpr std::uint64_t kSectionRows = 16;

struct Options {
    std::string data_dir;
    std::string stats_id;
    std::string run_id;
    std::string hp_check_id;
    bool toy = false;
    std::string direction;
    std::string truth;
    std::size_t events = 0;
    std::string acceptance_id;
    std::uint64_t seed = kProjectSeed;
    std::uint64_t split_seed = kProjectSeed;
    st::CaseBResolution resolution = kDefaultResolution;
    unsigned threads = kDefaultThreads;
    bool pilot = false;
};

// Exactly one source: a Geant4 run, or the toy with its direction, truth and size.
void requireOneSource(const Options& o) {
    if (o.toy) {
        if (!o.run_id.empty() || !o.hp_check_id.empty()) {
            throw apps::CliError("--toy excludes --run-id and --hp-from-check");
        }
        if (o.direction.empty() || o.truth.empty() || o.events == 0) {
            throw apps::CliError("--toy needs --direction, --truth and a positive --events");
        }
        static_cast<void>(apps::parseDirection(o.direction));
        if (apps::parseTruth(o.truth).mode != io::SpectrumMode::Flat) {
            throw apps::CliError("--truth: the Case B toy draws E_n from a flat spectrum");
        }
        return;
    }
    if (o.run_id.empty()) {
        throw apps::CliError("give either --run-id or --toy");
    }
    if (!o.direction.empty() || !o.truth.empty() || o.events != 0) {
        throw apps::CliError("--direction, --truth and --events go with --toy only");
    }
}

json directionJson(const Direction& d) {
    return json::array({d.x, d.y, d.z});
}

json domainJson(const st::EnergyDomain& d) {
    return {{"En_min_MeV", d.En_min}, {"En_max_MeV", d.En_max}, {"points", d.points}};
}

json truthJson(const io::SourceSpectrum& t) {
    json out = {{"mode", io::spectrumModeName(t.mode)}};
    if (t.mode == io::SpectrumMode::Mono) {
        out["En_MeV"] = t.En_MeV;
    } else {
        out["En_min_MeV"] = t.En_min_MeV;
        out["En_max_MeV"] = t.En_max_MeV;
    }
    return out;
}

json ensembleJson(const st::CaseBEnsemble& e) {
    return {{"entries", e.entries},
            {"error_mean_rad", e.error_mean},
            {"error_mean_error_rad", e.error_mean_error},
            {"error_rms_rad", e.error_rms},
            {"sigma_mean_rad", e.sigma_mean},
            {"sigma_mean_error_rad", e.sigma_mean_error},
            {"pull_rms", e.pull_rms},
            {"scaled_pull_median", e.scaled_pull_median},
            {"scaled_pull_tail_fraction", e.scaled_pull_tail_fraction},
            {"truth_in_cap_fraction", e.truth_in_cap_fraction},
            {"coverage", e.coverage}};
}

json biasJson(const st::DirectionBias& b) {
    return {{"basis", "e_theta, e_phi at the truth"},
            {"mean_rad", b.mean},
            {"error_rad", b.error},
            {"norm_rad", std::hypot(b.mean[0], b.mean[1])}};
}

class Stopwatch {
public:
    explicit Stopwatch(std::string name) : name_(std::move(name)) {}
    void report() const {
        const std::chrono::duration<double> elapsed = Clock::now() - start_;
        std::cout << "  " << name_ << ": " << elapsed.count() << " s\n" << std::flush;
    }

private:
    std::string name_;
    Clock::time_point start_ = Clock::now();
};

// The theta_p acceptance read from acceptance.json, and where it came from.
struct Acceptance {
    std::optional<st::AcceptanceRatio> ratio;
    json info; // null without --acceptance
};

Acceptance loadAcceptance(const fs::path& root, const std::string& id) {
    Acceptance out;
    if (id.empty()) {
        return out;
    }
    const json record = io::readJsonFile(riptide::data::statsDir(root, id) /
                                         std::string{io::kAcceptanceRecordName});
    const io::AcceptanceTable table = io::acceptanceTableFromRecord(record);
    out.ratio.emplace(st::AcceptanceHistogram{table.En_edges, table.theta_edges, table.counts});
    out.info = {{"stats_id", id},
                {"acceptance_version", record.at("acceptance_version")},
                {"run_id", record.at("source").at("run_id")},
                {"half", record.at("half")},
                {"entries", record.at("entries")},
                {"En_bins", table.En_edges.size() - 1},
                {"theta_bins", table.theta_edges.size() - 1},
                {"ratio_max", out.ratio->maximum()}};
    return out;
}

// The model and grids of 7a: E_n prior flat on the default domain, and with
// --acceptance the selection weights of D40 on the same grids.
struct Setup {
    st::EnergyDomain domain;
    st::Grid En = st::energyGrid(domain);
    st::Grid theta = st::recoilAngleGrid(st::kTrackAngleGridPoints);
    st::CaseBModel model;
    std::vector<Direction> coarse_grid = st::sphereGrid(st::kSphereGridPixels);

    Setup(const st::CaseBResolution& r, const st::AcceptanceRatio* acceptance)
        : model(En, theta, st::Grid{0.0, riptide::units::pi, st::kTrackObservedGridPoints},
                st::flatEnergyLogPrior(En), r,
                acceptance == nullptr ? std::vector<double>{}
                                      : st::caseBAcceptanceLogWeights(*acceptance, En, theta)) {}
};

// The pool of events, its true direction and spectrum, and what it came from.
struct Pool {
    st::CaseBPseudoSample sample;
    Direction omega_n;
    io::SourceSpectrum truth;
    json source;
};

Pool geant4Pool(const fs::path& root, const Options& o) {
    const std::optional<std::string> check =
        o.hp_check_id.empty() ? std::nullopt : std::optional<std::string>{o.hp_check_id};
    const io::RunInfo run = io::readRunInfo(root, o.run_id, check);
    const json run_record = io::readJsonFile(run.run_dir / io::kRunRecordName);
    const io::Selection selection = io::Selection::parse(kCaseBCuts);
    io::CutFlow flow(selection.cuts().size());
    Pool pool;
    pool.truth = io::sourceSpectrumFromRunRecord(run_record);
    pool.omega_n = apps::fixedDirectionFromRunRecord(run_record);
    {
        riptide::EventSample events;
        {
            // The full columns are needed only for the selection.
            const io::EventColumns columns = io::readEventColumns(run, 0);
            events = io::selectEvents(columns, selection, flow);
        }
        apps::requireFixedDirection(events, pool.omega_n);
        io::convertHpEnergies(events, run.hp.target_mass_MeV);
        io::applyResolution(events, io::Resolution{o.resolution.sigma_Ep, o.resolution.sigma_theta},
                            o.seed);
        pool.sample = apps::caseBSampleFromEvents(events);
    }
    pool.source = {{"kind", "geant4"},
                   {"run_id", run.run_id},
                   {"run_seed", run.seed},
                   {"run_record_version", run.run_record_version},
                   {"schema_version", run.schema_version},
                   {"run_git_describe", run.git_describe},
                   {"hp_mass",
                    {{"source", io::hpMassSourceName(run.hp.source)},
                     {"check_id", run.hp.check_id},
                     {"target_mass_MeV", run.hp.target_mass_MeV}}},
                   {"cuts", selection.toString()},
                   {"events_read", flow.total()},
                   {"events_selected", flow.selected()},
                   {"cut_flow", io::cutsJson(selection, flow)}};
    return pool;
}

Pool toyPool(const Setup& s, const st::AcceptanceRatio* acceptance, const Options& o) {
    Pool pool;
    pool.truth = apps::parseTruth(o.truth);
    pool.omega_n = apps::parseDirection(o.direction);
    const st::EnergyDomain truth{pool.truth.En_min_MeV, pool.truth.En_max_MeV, s.domain.points};
    pool.sample = acceptance == nullptr
                      ? st::caseBPseudoSample(truth, pool.omega_n, o.events, o.resolution, o.seed)
                      : st::caseBPseudoSample(truth, pool.omega_n, o.events, o.resolution, o.seed,
                                              *acceptance);
    pool.source = {{"kind", "toy"}, {"events", o.events}};
    return pool;
}

// The disjoint experiments of one row, or nothing when fewer than two fit:
// the row is then recorded as skipped, with what was planned and available.
struct RowSplit {
    std::size_t events = 0;
    std::size_t planned = 0;
    std::size_t available = 0; // events left for this row
    std::vector<std::vector<std::size_t>> experiments;
};

// Rows of the same N on one permutation of the pool, in order, so no event
// enters two of them; each takes what the previous ones left.
std::vector<RowSplit> splitRows(std::size_t pool, std::size_t events,
                                const std::vector<std::size_t>& planned, std::uint64_t section,
                                const Options& o) {
    std::vector<RowSplit> out;
    std::size_t left = pool;
    std::size_t total = 0;
    for (const std::size_t m_planned : planned) {
        RowSplit row;
        row.events = events;
        row.planned = apps::pilotScaled(m_planned, o.pilot);
        row.available = left;
        const std::size_t m = apps::rowExperiments(m_planned, left, events, o.pilot);
        row.experiments.resize(m);
        left -= m * events;
        total += m;
        out.push_back(std::move(row));
    }
    riptide::random::Engine engine = riptide::random::perEventEngine(o.split_seed, section);
    const std::vector<std::vector<std::size_t>> blocks =
        st::disjointExperiments(pool, events, total, engine);
    std::size_t next = 0;
    for (RowSplit& row : out) {
        for (std::vector<std::size_t>& experiment : row.experiments) {
            experiment = blocks.at(next++);
        }
    }
    return out;
}

struct RowResult {
    RowSplit split;
    std::optional<st::CaseBEnsemble> ensemble; // empty for a skipped row
    json record;
};

RowResult runRow(const Setup& s, const Pool& pool, RowSplit split, const Options& o,
                 const std::string& name) {
    const Stopwatch watch(name + " N = " + std::to_string(split.events) +
                          ", M = " + std::to_string(split.experiments.size()));
    RowResult out;
    out.record = {{"events", split.events},
                  {"experiments_planned", split.planned},
                  {"experiments", split.experiments.size()}};
    if (split.experiments.empty()) {
        out.record["experiments_available"] =
            st::maxDisjointExperiments(split.available, split.events);
        out.record["skipped"] = "too few events";
    } else {
        const std::vector<st::CaseBExperimentResult> results = st::runCaseBPoolRow(
            s.model, pool.sample, pool.omega_n, split.experiments, s.coarse_grid, o.threads);
        out.ensemble = st::summariseCaseBEnsemble(results);
        const double root_n = std::sqrt(static_cast<double>(split.events));
        out.record["skipped"] = nullptr;
        out.record["ensemble"] = ensembleJson(*out.ensemble);
        out.record["bias"] = biasJson(st::directionBias(results));
        out.record["error_rms_sqrt_N_rad"] = out.ensemble->error_rms * root_n;
    }
    out.split = std::move(split);
    watch.report();
    return out;
}

// log-log slope of the rms error over the analysed rows from kSlopeMinEvents
// up; null with fewer than two of them.
json slopeOf(const std::vector<const RowResult*>& rows) {
    std::vector<double> events;
    std::vector<double> rms;
    for (const RowResult* r : rows) {
        if (r->ensemble && r->split.events >= kSlopeMinEvents) {
            events.push_back(static_cast<double>(r->split.events));
            rms.push_back(r->ensemble->error_rms);
        }
    }
    if (events.size() < 2) {
        return nullptr;
    }
    return st::logLogSlope(events, rms);
}

int run(const Options& o) {
    namespace data = riptide::data;
    const std::optional<std::string> cli_root =
        o.data_dir.empty() ? std::nullopt : std::optional<std::string>{o.data_dir};
    const fs::path root = data::resolveDataRoot(cli_root);
    const fs::path out_dir = data::statsDir(root, o.stats_id);
    if (fs::exists(out_dir)) {
        throw data::DataPathError("stats directory already exists: " + out_dir.string());
    }
    data::requireWritableRoot(root, kOutputBytes, data::kFreeSpaceMargin, false);

    std::cout << kAppName << ' ' << riptide::kProjectVersion << " (" << riptide::kGitDescribe
              << "), " << o.threads << " threads" << (o.toy ? ", toy" : ", geant4")
              << (o.pilot ? ", pilot" : "") << '\n';
    const Acceptance acceptance = loadAcceptance(root, o.acceptance_id);
    const st::AcceptanceRatio* ratio = acceptance.ratio ? &*acceptance.ratio : nullptr;
    const Setup setup(o.resolution, ratio);
    const Pool pool = o.toy ? toyPool(setup, ratio, o) : geant4Pool(root, o);
    std::cout << "  pool: " << pool.sample.size() << " events\n";

    std::vector<RowResult> rows;
    for (std::size_t r = 0; r + 1 < kChecklistRows.size(); ++r) {
        const SizedRow row = kChecklistRows[r];
        std::vector<RowSplit> split =
            splitRows(pool.sample.size(), row.events, {row.experiments}, kSectionRows + r, o);
        rows.push_back(runRow(setup, pool, std::move(split.front()), o, "row"));
    }
    const SizedRow last = kChecklistRows.back();
    std::vector<RowSplit> last_split =
        splitRows(pool.sample.size(), last.events, {last.experiments, kRepeatExperiments},
                  kSectionRows + kChecklistRows.size() - 1, o);
    rows.push_back(runRow(setup, pool, std::move(last_split[0]), o, "row"));
    const RowResult repeat = runRow(setup, pool, std::move(last_split[1]), o, "repeat");

    std::vector<const RowResult*> slope_rows;
    json rows_json = json::array();
    for (const RowResult& r : rows) {
        slope_rows.push_back(&r);
        rows_json.push_back(r.record);
    }
    std::vector<const RowResult*> slope_rows_with_repeat(slope_rows.begin(), slope_rows.end() - 1);
    slope_rows_with_repeat.push_back(&repeat);
    // A skipped repeat would silently drop the largest N from the fit: null.
    const json slope_with_repeat =
        repeat.ensemble ? slopeOf(slope_rows_with_repeat) : json(nullptr);

    json record;
    record["case_b_version"] = kCaseBVersion;
    record["stats_id"] = o.stats_id;
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    record["source"] = pool.source;
    record["truth"] = {{"omega_n", directionJson(pool.omega_n)},
                       {"spectrum", truthJson(pool.truth)}};
    // The thread count is left out: it does not change a number (D44).
    record["settings"] = {{"seed", o.seed},
                          {"split_seed", o.split_seed},
                          {"sigma_Ep_MeV", o.resolution.sigma_Ep},
                          {"sigma_theta_rad", o.resolution.sigma_theta},
                          {"pilot", o.pilot},
                          {"levels", st::kDirectionCredibleLevels},
                          {"En_domain", domainJson(setup.domain)},
                          {"theta_points", st::kTrackAngleGridPoints},
                          {"theta_obs_points", st::kTrackObservedGridPoints},
                          {"sphere_pixels", st::kSphereGridPixels},
                          {"cap_pixels", st::kCapGridPixels},
                          {"cap_delta_log", st::kCapDeltaLog},
                          {"acceptance", acceptance.info}};
    record["rows"] = rows_json;
    record["repeat"] = repeat.record;
    record["contraction"] = {{"slope_min_events", kSlopeMinEvents},
                             {"slope", slopeOf(slope_rows)},
                             {"slope_with_repeat", slope_with_repeat}};

    fs::create_directories(out_dir);
    io::writeFileAtomically(out_dir / kCaseBName, record.dump(2) + "\n");
    std::cout << kAppName << ": output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    Options o;
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(o.stats_id, "id")["--stats-id"]("name of the output directory").required() |
        lyra::opt(o.data_dir, "dir")["--data-dir"]("data root (default: the environment)") |
        lyra::opt(o.run_id, "id")["--run-id"]("Geant4 source: the riptide_sim run") |
        lyra::opt(o.hp_check_id,
                  "id")["--hp-from-check"]("check with the HP mass (needed by run.json v1)") |
        lyra::opt(o.toy)["--toy"]("toy source, with --direction, --truth and --events") |
        lyra::opt(o.direction, "x,y,z")["--direction"]("toy neutron direction, a unit vector") |
        lyra::opt(o.truth, "spec")["--truth"]("toy spectrum flat:<min>:<max>, MeV") |
        lyra::opt(o.events, "n")["--events"]("events of the toy pool") |
        lyra::opt(o.seed, "n")["--seed"]("seed of the smearing or the toy (default 20260907)") |
        lyra::opt(o.split_seed, "n")["--split-seed"]("seed of the split (default 20260907)") |
        lyra::opt(o.resolution.sigma_Ep, "MeV")["--sigma-Ep"]("Ep resolution (default 0.10)") |
        lyra::opt(o.resolution.sigma_theta,
                  "rad")["--sigma-theta"]("track resolution per component (default 0.08)") |
        lyra::opt(o.acceptance_id,
                  "id")["--acceptance"]("stats id of acceptance.json: selection weights (D40)") |
        lyra::opt(o.threads, "n")["--threads"]("threads across experiments (default 4)") |
        lyra::opt(o.pilot)["--pilot"]("counts / 10, same N: to measure time and memory");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    try {
        if (!result) {
            throw apps::CliError(result.message());
        }
        if (o.threads == 0) {
            throw apps::CliError("--threads must be at least 1");
        }
        requireOneSource(o);
    } catch (const apps::CliError& e) {
        std::cerr << kAppName << ": " << e.what() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    try {
        return run(o);
    } catch (const std::exception& e) {
        std::cerr << kAppName << ": " << e.what() << '\n';
        return kExitRuntime;
    }
}
