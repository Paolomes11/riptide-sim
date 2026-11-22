#include "riptide/sim/run_reader.hpp"

#include <G4RootAnalysisReader.hh>

#include <array>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

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

// One event file as a RowSource. The first n_columns schema columns (those of
// the file's schema version) are bound to buffers owned here, so the source
// must not move once bound; the later columns of the row stay 0 of their type.
class FileRowSource final : public RowSource {
public:
    FileRowSource(std::filesystem::path file, std::size_t n_columns)
        : file_(std::move(file)), n_columns_(n_columns) {
        if (!std::filesystem::is_regular_file(file_)) {
            throw RunReadError("missing event file " + file_.string());
        }
        auto* reader = G4RootAnalysisReader::Instance();
        id_ = reader->GetNtuple(G4String(schema::kEventTreeName), file_.string());
        if (id_ == kInvalidNtupleId) {
            throw RunReadError("cannot read ntuple from " + file_.string());
        }
        for (std::size_t i = 0; i < n_columns_; ++i) {
            const auto& spec = columnSpec(i);
            const G4String col(spec.name);
            const bool bound = spec.type == ColumnType::Int32
                                   ? reader->SetNtupleIColumn(id_, col, ints_.at(i))
                                   : reader->SetNtupleDColumn(id_, col, doubles_.at(i));
            if (!bound) {
                throw RunReadError("cannot bind column " + col + " in " + file_.string());
            }
        }
    }

    bool next(EventRow& row) override {
        if (!G4RootAnalysisReader::Instance()->GetNtupleRow(id_)) {
            return false;
        }
        for (std::size_t i = 0; i < kEventColumns; ++i) {
            const bool int_column = columnSpec(i).type == ColumnType::Int32;
            if (i < n_columns_) {
                row.at(i) = int_column ? ColumnValue{ints_.at(i)} : ColumnValue{doubles_.at(i)};
            } else {
                row.at(i) = int_column ? ColumnValue{std::int32_t{0}} : ColumnValue{0.0};
            }
        }
        return true;
    }

    [[nodiscard]] std::string name() const override {
        return file_.string();
    }

private:
    std::filesystem::path file_;
    std::size_t n_columns_;
    G4int id_ = kInvalidNtupleId;
    std::array<G4int, kEventColumns> ints_{};
    std::array<G4double, kEventColumns> doubles_{};
};

// The head of one source during the merge: its current row, if any.
struct MergeHead {
    RowSource* source = nullptr;
    EventRow row;
    std::int32_t event_id = 0;
    bool live = false;
};

void advance(MergeHead& head, std::size_t event_id_col) {
    const bool had_row = head.live;
    const std::int32_t previous = head.event_id;
    head.live = head.source->next(head.row);
    if (!head.live) {
        return;
    }
    head.event_id = std::get<std::int32_t>(head.row.at(event_id_col));
    if (had_row && head.event_id <= previous) {
        throw RunReadError("event_id " + std::to_string(head.event_id) + " after " +
                           std::to_string(previous) + " in " + head.source->name() +
                           ": rows are not in ascending event_id");
    }
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
        const auto record_version = record.at("run_record_version").get<int>();
        if (record_version < kOldestReadableRunRecordVersion ||
            record_version > kRunRecordVersion) {
            throw RunReadError("unknown run_record_version in " + path.string());
        }
        const auto version = record.at("schema_version").get<std::int32_t>();
        if (version < schema::kOldestReadableSchemaVersion || version > schema::kSchemaVersion) {
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

std::size_t mergeByEventId(const std::vector<RowSource*>& sources,
                           const std::function<void(const EventRow&)>& visit) {
    const std::size_t event_id_col = columnIndex("event_id");
    std::vector<MergeHead> heads(sources.size());
    for (std::size_t i = 0; i < sources.size(); ++i) {
        heads[i].source = sources[i];
        advance(heads[i], event_id_col);
    }
    std::size_t count = 0;
    // A linear scan for the smallest head: there is one source per worker
    // thread, a handful, so a heap would not pay for itself.
    for (;;) {
        MergeHead* lowest = nullptr;
        for (auto& head : heads) {
            if (head.live && (lowest == nullptr || head.event_id < lowest->event_id)) {
                lowest = &head;
            }
        }
        if (lowest == nullptr) {
            return count;
        }
        visit(lowest->row);
        ++count;
        advance(*lowest, event_id_col);
    }
}

std::size_t forEachEventRow(const std::filesystem::path& run_dir, const json& record,
                            const std::function<void(const EventRow&)>& visit) {
    const std::size_t n_columns = columnCount(record.at("schema_version").get<std::int32_t>());
    const ReaderFilesGuard guard;
    // Every file stays open for the whole merge; the Geant4 reader keeps one
    // ntuple description, with its own row cursor, per GetNtuple call.
    std::vector<std::unique_ptr<FileRowSource>> files;
    std::vector<RowSource*> sources;
    for (const auto& name : record.at("event_files")) {
        files.push_back(
            std::make_unique<FileRowSource>(run_dir / name.get<std::string>(), n_columns));
        sources.push_back(files.back().get());
    }
    return mergeByEventId(sources, visit);
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
