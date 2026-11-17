#include "riptide/sim/ntuple_writer.hpp"

#include <G4AnalysisManager.hh>
#include <G4Exception.hh>
#include <G4String.hh>

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <variant>

namespace riptide::sim {

namespace {

// Geant4 prints nothing about bookings and file operations at this level.
constexpr G4int kAnalysisVerbosity = 0;

void requireOk(bool ok, const char* what, const std::string& detail) {
    if (!ok) {
        // A worker cannot propagate a C++ exception to main; G4Exception with
        // FatalException stops the whole run with a message instead.
        G4Exception("riptide::sim::NtupleWriter", "RIPTIDE_NTUPLE", FatalException,
                    (std::string(what) + ": " + detail).c_str());
    }
}

} // namespace

NtupleWriter::NtupleWriter() {
    auto* manager = G4AnalysisManager::Instance();
    manager->SetVerboseLevel(kAnalysisVerbosity);
    manager->SetNtupleMerging(false);
    const std::string tree_name(schema::kEventTreeName);
    ntuple_id_ = manager->CreateNtuple(tree_name, tree_name);
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        const auto& spec = columnSpec(i);
        const G4String name(std::string(spec.name));
        column_ids_.at(i) = spec.type == schema::ColumnType::Int32
                                ? manager->CreateNtupleIColumn(ntuple_id_, name)
                                : manager->CreateNtupleDColumn(ntuple_id_, name);
    }
    manager->FinishNtuple(ntuple_id_);
}

void NtupleWriter::open(const std::filesystem::path& path) {
    requireOk(G4AnalysisManager::Instance()->OpenFile(path.string()), "cannot open", path.string());
}

void NtupleWriter::fill(const EventRow& row) {
    auto* manager = G4AnalysisManager::Instance();
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        const G4int column = column_ids_.at(i);
        const bool ok = std::visit(
            [&](auto value) {
                if constexpr (std::is_same_v<decltype(value), std::int32_t>) {
                    return manager->FillNtupleIColumn(ntuple_id_, column, value);
                } else {
                    return manager->FillNtupleDColumn(ntuple_id_, column, value);
                }
            },
            row.at(i));
        requireOk(ok, "cannot fill column", std::string(columnSpec(i).name));
    }
    requireOk(manager->AddNtupleRow(ntuple_id_), "cannot add row", "events ntuple");
}

void NtupleWriter::close() {
    auto* manager = G4AnalysisManager::Instance();
    requireOk(manager->Write(), "cannot write", "events ntuple");
    requireOk(manager->CloseFile(), "cannot close", "events ntuple");
}

} // namespace riptide::sim
