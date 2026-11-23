#include "riptide/io/analysis_record.hpp"

#include <fstream>
#include <sstream>
#include <system_error>

#include "riptide/common/git_info.hpp"
#include "riptide/common/version.hpp"

namespace riptide::io {

namespace {

constexpr std::string_view kTemporarySuffix = ".tmp";
constexpr int kRoundTripDigits = 17;

} // namespace

nlohmann::json cutsJson(const Selection& selection, const CutFlow& cut_flow) {
    if (cut_flow.cuts() != selection.cuts().size()) {
        throw IoError("the cut flow does not belong to the selection");
    }
    auto cuts = nlohmann::json::array();
    for (std::size_t k = 0; k < selection.cuts().size(); ++k) {
        const auto& spec = selection.cuts()[k];
        nlohmann::json cut;
        cut["name"] = std::string{cutName(spec.cut)};
        if (spec.cut == Cut::EvisMin) {
            cut["threshold_MeV"] = spec.threshold;
        }
        cut["passed"] = cut_flow.passed(k);
        cuts.push_back(cut);
    }
    return cuts;
}

nlohmann::json analysisRecordJson(const RunInfo& run, const AnalysisSettings& settings,
                                  const CutFlow& cut_flow, const nlohmann::json& summaries) {
    if (cut_flow.cuts() != settings.selection.cuts().size()) {
        throw IoError("the cut flow does not belong to the selection");
    }
    nlohmann::json doc;
    doc["analysis_record_version"] = kAnalysisRecordVersion;
    doc["analysis_id"] = settings.analysis_id;
    doc["project_version"] = std::string{kProjectVersion};
    doc["git_describe"] = std::string{kGitDescribe};

    nlohmann::json input;
    input["run_id"] = run.run_id;
    input["run_record_version"] = run.run_record_version;
    input["schema_version"] = run.schema_version;
    input["sim_git_describe"] = run.git_describe;
    input["events_in_run"] = run.events;
    input["max_events"] = settings.max_events;
    doc["input"] = input;

    nlohmann::json hp;
    hp["source"] = std::string{hpMassSourceName(run.hp.source)};
    if (run.hp.source == HpMassSource::Check) {
        hp["check_id"] = run.hp.check_id;
    }
    hp["file"] = run.hp.file;
    hp["sha256"] = run.hp.sha256;
    hp["awr"] = run.hp.awr;
    hp["target_mass_MeV"] = run.hp.target_mass_MeV;
    doc["hp_target_mass"] = hp;

    nlohmann::json response;
    response["seed"] = settings.seed;
    response["sigma_Ep_MeV"] = settings.resolution.sigma_Ep;
    response["sigma_theta_rad"] = settings.resolution.sigma_theta;
    response["per_event_engine"] = "mt19937_64 seeded by splitmix64(seed, event_id)";
    doc["response"] = response;

    nlohmann::json selection;
    selection["text"] = settings.selection.toString();
    selection["events_read"] = cut_flow.total();
    selection["events_selected"] = cut_flow.selected();
    selection["cuts"] = cutsJson(settings.selection, cut_flow);
    doc["selection"] = selection;

    doc["summaries"] = summaries;
    return doc;
}

std::string cutFlowCsv(const Selection& selection, const CutFlow& cut_flow) {
    if (cut_flow.cuts() != selection.cuts().size()) {
        throw IoError("the cut flow does not belong to the selection");
    }
    std::ostringstream out;
    out.precision(kRoundTripDigits);
    const auto total = static_cast<double>(cut_flow.total());
    const auto fraction = [&](std::uint64_t count) {
        return cut_flow.total() == 0 ? 0.0 : static_cast<double>(count) / total;
    };
    out << "index,cut,threshold_MeV,passed,fraction\n";
    out << "0,all,," << cut_flow.total() << ',' << fraction(cut_flow.total()) << '\n';
    for (std::size_t k = 0; k < selection.cuts().size(); ++k) {
        const auto& spec = selection.cuts()[k];
        out << k + 1 << ',' << cutName(spec.cut) << ',';
        if (spec.cut == Cut::EvisMin) {
            out << spec.threshold;
        }
        out << ',' << cut_flow.passed(k) << ',' << fraction(cut_flow.passed(k)) << '\n';
    }
    return out.str();
}

void writeFileAtomically(const std::filesystem::path& target, const std::string& content) {
    if (std::filesystem::exists(target)) {
        throw IoError("file already exists: " + target.string());
    }
    auto temporary = target;
    temporary += kTemporarySuffix;
    {
        std::ofstream out(temporary, std::ios::out | std::ios::trunc);
        out << content;
        out.close();
        if (!out) {
            throw IoError("cannot write " + temporary.string());
        }
    }
    // rename is atomic within one filesystem.
    std::error_code error;
    std::filesystem::rename(temporary, target, error);
    if (error) {
        throw IoError("cannot rename " + temporary.string() + " to " + target.string() + ": " +
                      error.message());
    }
}

} // namespace riptide::io
