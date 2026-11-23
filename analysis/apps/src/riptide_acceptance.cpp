// theta_p acceptance of Case A (D40): reads one riptide_sim run, applies the
// Case A cut chain and histograms truth_En and truth_theta_p of the selected
// events, 1 MeV by 5 degrees. It writes acceptance.json under the stats
// directory of the data root, for riptide_case_a_geant4 --acceptance. The
// table should come from a run independent of the ones it is used on, so its
// fluctuations do not enter twice. The file holds no wall time: the same run
// and settings give the same bytes.

#include <lyra/lyra.hpp>

#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_sample.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/units.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/acceptance_record.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace {

namespace fs = std::filesystem;
namespace apps = riptide::apps;
namespace io = riptide::io;

constexpr const char* kAppName = "riptide_acceptance";
constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;
// Generous upper bound on the output: one small JSON file.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB
constexpr double kHalf = 0.5;

// The cut chain of Case A (D30, D34), the same as riptide_case_a_geant4.
constexpr const char* kCaseACuts = "single_np,has_proton,contained";

struct Options {
    std::string data_dir;
    std::string stats_id;
    std::string run_id;
    std::string half = "all";
};

// E_n edges over the source range. A mono run has one bin centred on its
// energy, so the bin centre, where the ratio is evaluated, is the energy
// itself; its width only places the events.
std::vector<double> energyEdges(const io::SourceSpectrum& source) {
    if (source.mode == io::SpectrumMode::Mono) {
        const double half_step = kHalf * io::kAcceptanceEnStep;
        return {source.En_MeV - half_step, source.En_MeV + half_step};
    }
    return io::steppedEdges(source.En_min_MeV, source.En_max_MeV, io::kAcceptanceEnStep);
}

int run(const Options& o, io::EventHalf half) {
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
              << ")\n";
    // truth_En and truth_theta_p are read as written: neither depends on the
    // HP target-mass convention, so no HP mass is needed.
    const io::RunInfo run = io::readRunInfo(root, o.run_id, std::nullopt);
    const io::SourceSpectrum source =
        io::sourceSpectrumFromRunRecord(io::readJsonFile(run.run_dir / io::kRunRecordName));
    const io::Selection selection = io::Selection::parse(kCaseACuts);
    io::CutFlow flow(selection.cuts().size());
    riptide::EventSample events;
    {
        const io::EventColumns columns = io::readEventColumns(run, 0);
        events = io::selectEvents(columns, selection, flow);
    }
    const io::AcceptanceTable table = io::fillAcceptanceTable(
        events, energyEdges(source),
        io::uniformEdges(0.0, kHalf * riptide::units::pi, io::kAcceptanceThetaBins), half);

    nlohmann::json record = io::acceptanceRecordJson(o.stats_id, run, selection, flow, half, table);
    record["project_version"] = riptide::kProjectVersion;
    record["git_describe"] = riptide::kGitDescribe;
    fs::create_directories(out_dir);
    io::writeFileAtomically(out_dir / io::kAcceptanceRecordName, record.dump(2) + "\n");
    std::cout << "  selected " << flow.selected() << " of " << flow.total() << ", entries "
              << record.at("entries").get<double>() << '\n'
              << kAppName << ": output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    bool show_help = false;
    Options o;
    const auto cli =
        lyra::cli() | lyra::help(show_help) |
        lyra::opt(o.stats_id, "id")["--stats-id"]("name of the output directory").required() |
        lyra::opt(o.run_id, "id")["--run-id"]("the riptide_sim run").required() |
        lyra::opt(o.data_dir, "dir")["--data-dir"]("data root (default: the environment)") |
        lyra::opt(o.half, "all|even|odd")["--half"]("events by event_id parity (default all)");
    const auto result = cli.parse(lyra::args(argc, argv));
    if (show_help) {
        std::cout << cli << '\n';
        return kExitOk;
    }
    io::EventHalf half = io::EventHalf::All;
    try {
        if (!result) {
            throw apps::CliError(result.message());
        }
        half = io::parseEventHalf(o.half);
    } catch (const std::exception& e) {
        std::cerr << kAppName << ": " << e.what() << "\n\n" << cli << '\n';
        return kExitUsage;
    }
    try {
        return run(o, half);
    } catch (const std::exception& e) {
        std::cerr << kAppName << ": " << e.what() << '\n';
        return kExitRuntime;
    }
}
