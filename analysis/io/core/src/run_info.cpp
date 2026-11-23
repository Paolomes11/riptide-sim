#include "riptide/io/run_info.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <utility>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"

namespace riptide::io {

namespace {

using nlohmann::json;

const json& member(const json& doc, const char* key, const std::string& where) {
    if (!doc.is_object() || !doc.contains(key)) {
        throw IoError(where + ": missing key '" + key + "'");
    }
    return doc.at(key);
}

std::string stringMember(const json& doc, const char* key, const std::string& where) {
    const auto& value = member(doc, key, where);
    if (!value.is_string()) {
        throw IoError(where + ": '" + key + "' is not a string");
    }
    return value.get<std::string>();
}

std::uint64_t unsignedMember(const json& doc, const char* key, const std::string& where) {
    const auto& value = member(doc, key, where);
    if (!value.is_number_unsigned() &&
        !(value.is_number_integer() && value.get<std::int64_t>() >= 0)) {
        throw IoError(where + ": '" + key + "' is not a non-negative integer");
    }
    return value.get<std::uint64_t>();
}

int intMember(const json& doc, const char* key, const std::string& where) {
    const auto& value = member(doc, key, where);
    if (!value.is_number_integer()) {
        throw IoError(where + ": '" + key + "' is not an integer");
    }
    return value.get<int>();
}

double positiveMember(const json& doc, const char* key, const std::string& where) {
    const auto& value = member(doc, key, where);
    if (!value.is_number()) {
        throw IoError(where + ": '" + key + "' is not a number");
    }
    const auto number = value.get<double>();
    if (!std::isfinite(number) || number <= 0.0) {
        throw IoError(where + ": '" + key + "' is not positive and finite");
    }
    return number;
}

// An event file name must stay inside the run directory.
bool isPlainFileName(const std::string& name) {
    return data::isValidId(name);
}

} // namespace

std::string_view hpMassSourceName(HpMassSource source) noexcept {
    switch (source) {
    case HpMassSource::RunRecord:
        return "run.json";
    case HpMassSource::Check:
        return "checks.json";
    }
    return "unknown";
}

HpMass hpMassFromBlock(const json& block, HpMassSource source, std::string check_id) {
    const std::string where = "hp_elastic_H1";
    HpMass hp;
    hp.awr = positiveMember(block, "awr", where);
    hp.target_mass_MeV = positiveMember(block, "target_mass_MeV", where);
    hp.file = stringMember(block, "file", where);
    hp.sha256 = stringMember(block, "sha256", where);
    hp.source = source;
    hp.check_id = std::move(check_id);
    const double expected = hp.awr * kinematics::kNeutronMass;
    if (std::abs(hp.target_mass_MeV - expected) > kHpMassRelTolerance * expected) {
        throw IoError(where + ": target_mass_MeV does not equal awr times the neutron mass");
    }
    return hp;
}

std::string_view spectrumModeName(SpectrumMode mode) noexcept {
    switch (mode) {
    case SpectrumMode::Mono:
        return "mono";
    case SpectrumMode::Flat:
        return "flat";
    }
    return "unknown";
}

SourceSpectrum sourceSpectrumFromRunRecord(const json& record) {
    const std::string record_where{kRunRecordName};
    const std::string config_where = record_where + ": config";
    const std::string where = config_where + ".source";
    const auto& config = member(record, "config", record_where);
    const auto& source = member(config, "source", config_where);
    const auto mode = stringMember(source, "mode", where);
    SourceSpectrum spectrum;
    if (mode == spectrumModeName(SpectrumMode::Mono)) {
        spectrum.mode = SpectrumMode::Mono;
        spectrum.En_MeV = positiveMember(source, "En_MeV", where);
    } else if (mode == spectrumModeName(SpectrumMode::Flat)) {
        spectrum.mode = SpectrumMode::Flat;
        spectrum.En_min_MeV = positiveMember(source, "En_min_MeV", where);
        spectrum.En_max_MeV = positiveMember(source, "En_max_MeV", where);
        if (spectrum.En_min_MeV >= spectrum.En_max_MeV) {
            throw IoError(where + ": En_min_MeV is not below En_max_MeV");
        }
    } else {
        throw IoError(where + ": unknown mode '" + mode + "'");
    }
    return spectrum;
}

RunInfo runInfoFromJson(const std::filesystem::path& run_dir, const json& record,
                        const std::optional<std::string>& hp_check_id,
                        const std::optional<json>& check) {
    const std::string where = (run_dir / kRunRecordName).string();
    RunInfo info;
    info.run_dir = run_dir;
    info.run_record_version = intMember(record, "run_record_version", where);
    if (info.run_record_version < kOldestRunRecordVersion ||
        info.run_record_version > kNewestRunRecordVersion) {
        throw IoError(where + ": unknown run_record_version " +
                      std::to_string(info.run_record_version));
    }
    info.run_id = stringMember(record, "run_id", where);
    if (info.run_id != run_dir.filename().string()) {
        throw IoError(where + ": run_id '" + info.run_id + "' does not match its directory");
    }
    info.schema_version = intMember(record, "schema_version", where);
    if (info.schema_version < schema::kOldestReadableSchemaVersion ||
        info.schema_version > schema::kSchemaVersion) {
        throw IoError(where + ": unreadable schema_version " + std::to_string(info.schema_version));
    }
    info.events = unsignedMember(record, "events", where);
    info.seed = unsignedMember(record, "seed", where);
    info.git_describe = stringMember(record, "git_describe", where);
    const auto& files = member(record, "event_files", where);
    if (!files.is_array() || files.empty()) {
        throw IoError(where + ": event_files is not a non-empty list");
    }
    for (const auto& name : files) {
        if (!name.is_string() || !isPlainFileName(name.get<std::string>())) {
            throw IoError(where + ": event_files holds an invalid file name");
        }
        info.event_files.push_back(run_dir / name.get<std::string>());
    }

    if (check.has_value() != hp_check_id.has_value()) {
        throw IoError("internal: a check document needs its check id and vice versa");
    }
    std::optional<HpMass> from_check;
    if (check) {
        const std::string check_where =
            "check '" + *hp_check_id + "' (" + std::string{kChecksName} + ")";
        if (stringMember(*check, "run_id", check_where) != info.run_id) {
            throw IoError(check_where + ": belongs to another run");
        }
        if (stringMember(*check, "check_id", check_where) != *hp_check_id) {
            throw IoError(check_where + ": check_id does not match its directory");
        }
        from_check = hpMassFromBlock(member(*check, "hp_elastic_H1", check_where),
                                     HpMassSource::Check, *hp_check_id);
    }

    if (info.run_record_version >= kFirstRunRecordWithHp) {
        info.hp =
            hpMassFromBlock(member(record, "hp_elastic_H1", where), HpMassSource::RunRecord, {});
        if (from_check) {
            // Both blocks are written by the same code from the same G4NDL file:
            // any difference means the run and the check saw different nuclear data.
            const bool same = record.at("hp_elastic_H1") == check->at("hp_elastic_H1");
            if (!same) {
                throw IoError("hp_elastic_H1 of run.json and of check '" + *hp_check_id +
                              "' differ");
            }
        }
    } else if (from_check) {
        info.hp = *from_check;
    } else {
        throw IoError(where + ": run_record_version " + std::to_string(info.run_record_version) +
                      " has no HP target mass; name a check of this run that reports it");
    }
    return info;
}

json readJsonFile(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw IoError("cannot open " + path.string());
    }
    try {
        return json::parse(in);
    } catch (const json::parse_error& error) {
        throw IoError("malformed JSON in " + path.string() + ": " + error.what());
    }
}

RunInfo readRunInfo(const std::filesystem::path& data_root, std::string_view run_id,
                    const std::optional<std::string>& hp_check_id) {
    const auto run_dir = data::simRunDir(data_root, run_id);
    const auto record = readJsonFile(run_dir / kRunRecordName);
    std::optional<json> check;
    if (hp_check_id) {
        check = readJsonFile(data::validationDir(data_root, run_id, *hp_check_id) / kChecksName);
    }
    auto info = runInfoFromJson(run_dir, record, hp_check_id, check);
    for (const auto& file : info.event_files) {
        if (!std::filesystem::is_regular_file(file)) {
            throw IoError("event file listed in run.json does not exist: " + file.string());
        }
    }
    return info;
}

} // namespace riptide::io
