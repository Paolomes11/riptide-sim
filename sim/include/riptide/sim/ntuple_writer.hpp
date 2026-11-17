#pragma once

#include <G4Types.hh>

#include <array>
#include <filesystem>

#include "riptide/sim/event_row.hpp"

namespace riptide::sim {

// Writes the event ntuple of one worker thread through the thread's Geant4
// analysis manager. Ntuple merging is off, so every worker writes its own file:
// Geant4 turns the path "<dir>/events.root" into "<dir>/events_t<thread id>.root".
// The columns are booked from event_schema.hpp, which stays the only source of
// names and types. Failures to open or write the file are fatal: a run whose
// output is incomplete must not look like a successful one.
class NtupleWriter {
public:
    // Books the ntuple; call once per worker thread, before the first open().
    NtupleWriter();

    void open(const std::filesystem::path& path);
    void fill(const EventRow& row);
    void close();

private:
    G4int ntuple_id_ = -1;
    std::array<G4int, kEventColumns> column_ids_{};
};

} // namespace riptide::sim
