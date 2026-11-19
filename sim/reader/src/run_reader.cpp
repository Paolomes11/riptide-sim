#include "riptide/sim/run_reader.hpp"

#include <G4RootAnalysisReader.hh>

#include <array>
#include <fstream>
#include <string>
#include <variant>

#include "riptide/common/event_schema.hpp"
#include "riptide/sim/run_record.hpp"

namespace riptide::sim {

namespace {

using nlohmann::json;
using schema::ColumnType;

// G4AnalysisReader's "no ntuple" id (kInvalidId in G4AnalysisUtilities.hh).
constexpr G4int kInvalidNtupleId = -1;

// Closes the reader's files however the reading ends, so a failed read does
// not leave a file open for the next run read in the same process.
class ReaderFilesGuard {
public:
    ReaderFilesGuard() = default;
    ReaderFilesGuard(const ReaderFilesGuard&) = delete;
    ReaderFilesGuard& operator=(const ReaderFilesGuard&) = delete;
    ReaderFilesGuard(ReaderFilesGuard&&) = delete;
    ReaderFilesGuard& operator=(ReaderFilesGuard&&) = delete;
    ~ReaderFilesGuard() {
        G4RootAnalysisReader::Instance()->CloseFiles();
    }
};

std::size_t forEachRowInFile(const std::filesystem::path& file,
                             const std::function<void(const EventRow&)>& visit) {
    if (!std::filesystem::is_regular_file(file)) {
        throw RunReadError("missing event file " + file.string());
    }
    auto* reader = G4RootAnalysisReader::Instance();
    const G4int id = reader->GetNtuple(G4String(schema::kEventTreeName), file.string());
    if (id == kInvalidNtupleId) {
        throw RunReadError("cannot read ntuple from " + file.string());
    }
    std::array<G4int, kEventColumns> ints{};
    std::array<G4double, kEventColumns> doubles{};
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        const auto& spec = columnSpec(i);
        const G4String name(spec.name);
        const bool bound = spec.type == ColumnType::Int32
                               ? reader->SetNtupleIColumn(id, name, ints.at(i))
                               : reader->SetNtupleDColumn(id, name, doubles.at(i));
        if (!bound) {
            throw RunReadError("cannot bind column " + name + " in " + file.string());
        }
    }
    std::size_t count = 0;
    EventRow row;
    while (reader->GetNtupleRow(id)) {
        for (std::size_t i = 0; i < kEventColumns; ++i) {
            row.at(i) = columnSpec(i).type == ColumnType::Int32 ? ColumnValue{ints.at(i)}
                                                                : ColumnValue{doubles.at(i)};
        }
        visit(row);
        ++count;
    }
    return count;
}

} // namespace

json readRunRecord(const std::filesystem::path& run_dir) {
    const auto path = run_dir / kRunRecordFileName;
    std::ifstream in(path);
    if (!in) {
        throw RunReadError("no " + std::string{kRunRecordFileName} + " in " + run_dir.string() +
                           " (interrupted or missing run)");
    }
    json record;
    try {
        record = json::parse(in);
        if (record.at("run_record_version").get<int>() != kRunRecordVersion) {
            throw RunReadError("unknown run_record_version in " + path.string());
        }
        if (record.at("schema_version").get<std::int32_t>() != schema::kSchemaVersion) {
            throw RunReadError("unknown schema_version in " + path.string());
        }
        const auto& files = record.at("event_files");
        if (!files.is_array()) {
            throw RunReadError("event_files is not a list in " + path.string());
        }
        for (const auto& name : files) {
            if (!name.is_string()) {
                throw RunReadError("event_files holds a non-string in " + path.string());
            }
        }
    } catch (const json::exception& e) {
        throw RunReadError("malformed " + path.string() + ": " + e.what());
    }
    return record;
}

std::size_t forEachEventRow(const std::filesystem::path& run_dir, const json& record,
                            const std::function<void(const EventRow&)>& visit) {
    const ReaderFilesGuard guard;
    std::size_t count = 0;
    for (const auto& name : record.at("event_files")) {
        count += forEachRowInFile(run_dir / name.get<std::string>(), visit);
    }
    return count;
}

RunData loadRun(const std::filesystem::path& run_dir) {
    RunData run;
    run.record = readRunRecord(run_dir);
    const std::size_t event_id_col = columnIndex("event_id");
    run.row_count = forEachEventRow(run_dir, run.record, [&run, event_id_col](const EventRow& row) {
        run.rows.emplace(std::get<std::int32_t>(row.at(event_id_col)), row);
    });
    return run;
}

std::int32_t intAt(const EventRow& row, std::string_view name) {
    return std::get<std::int32_t>(row.at(columnIndex(name)));
}

double doubleAt(const EventRow& row, std::string_view name) {
    return std::get<double>(row.at(columnIndex(name)));
}

} // namespace riptide::sim
