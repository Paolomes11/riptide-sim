// Prepares the sample of an analysis from one riptide_sim run: reads the
// event files, applies the cut chain, converts the HP proton energies to the
// physical ones and draws the detector resolution with a recorded seed. It
// writes only the record of what it did (analysis.json, cutflow.csv) and the
// residual summaries that show the conversion and the smearing behaved; the
// sample itself stays in memory, where the statistics stage will use it.

#include <nlohmann/json.hpp>

#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "riptide/apps/prepare_cli.hpp"
#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"
#include "riptide/io/analysis_record.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/response_summary.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace {

namespace fs = std::filesystem;
namespace apps = riptide::apps;
namespace io = riptide::io;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 2;
constexpr int kExitRuntime = 1;

// Generous upper bound on the output: two small text files whose size does
// not grow with the run.
constexpr std::uintmax_t kOutputBytes = std::uintmax_t{1} << 20U; // 1 MiB

void printVersion() {
    std::cout << "riptide_analysis_prepare " << riptide::kProjectVersion << " ("
              << riptide::kGitDescribe << ")\n"
              << "event schema " << riptide::schema::kSchemaVersion << '\n'
              << "analysis record version " << io::kAnalysisRecordVersion << '\n';
}

int prepare(const apps::PrepareOptions& opts) {
    namespace data = riptide::data;
    const fs::path root = data::resolveDataRoot(opts.data_dir);
    const fs::path out_dir = data::analysisDir(root, opts.run_id, opts.analysis_id);
    if (fs::exists(out_dir)) {
        throw data::DataPathError("analysis directory already exists: " + out_dir.string());
    }
    data::requireWritableRoot(root, kOutputBytes, data::kFreeSpaceMargin, false);

    const io::RunInfo run = io::readRunInfo(root, opts.run_id, opts.hp_check_id);
    io::CutFlow flow(opts.selection.cuts().size());
    riptide::EventSample sample;
    {
        // The full columns are needed only for the selection; free them before
        // the sample grows its reco columns.
        const io::EventColumns columns = io::readEventColumns(run, opts.max_events);
        sample = io::selectEvents(columns, opts.selection, flow);
    }
    io::convertHpEnergies(sample, run.hp.target_mass_MeV);
    io::applyResolution(sample, opts.resolution, opts.seed);

    nlohmann::json summaries;
    summaries["hp_conversion"] = io::hpConversionSummary(sample);
    summaries["resolution"] = io::resolutionSummary(sample, opts.resolution);

    io::AnalysisSettings settings;
    settings.analysis_id = opts.analysis_id;
    settings.seed = opts.seed;
    settings.resolution = opts.resolution;
    settings.selection = opts.selection;
    settings.max_events = opts.max_events;

    fs::create_directories(out_dir);
    io::writeFileAtomically(out_dir / io::kCutFlowName, io::cutFlowCsv(opts.selection, flow));
    // Last: its presence marks a complete analysis.
    io::writeFileAtomically(out_dir / io::kAnalysisRecordName,
                            io::analysisRecordJson(run, settings, flow, summaries).dump(2) + "\n");

    std::cout << "riptide_analysis_prepare: " << flow.total() << " events read, " << flow.selected()
              << " selected, output in " << out_dir << '\n';
    return kExitOk;
}

} // namespace

int main(int argc, char** argv) {
    apps::PrepareOptions opts;
    try {
        opts = apps::parsePrepareCommandLine(std::vector<std::string>(argv + 1, argv + argc));
    } catch (const apps::CliError& e) {
        std::cerr << "riptide_analysis_prepare: " << e.what() << "\n\n" << apps::prepareHelpText();
        return kExitUsage;
    }
    if (opts.show_help) {
        std::cout << apps::prepareHelpText();
        return kExitOk;
    }
    if (opts.show_version) {
        printVersion();
        return kExitOk;
    }
    try {
        return prepare(opts);
    } catch (const std::exception& e) {
        std::cerr << "riptide_analysis_prepare: " << e.what() << '\n';
        return kExitRuntime;
    }
}
