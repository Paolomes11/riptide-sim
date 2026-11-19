#pragma once

#include <G4UserRunAction.hh>

#include <filesystem>
#include <utility>

#include "riptide/sim/deposit_record.hpp"
#include "riptide/sim/ntuple_writer.hpp"
#include "riptide/sim/truth_record.hpp"

class G4Run;

namespace riptide::sim {

// Worker-thread run action. It owns the thread's truth and deposit recorders
// (shared by reference with the generator, stepping and event actions), the run summary
// printed at the end of the run and the thread's ntuple file in run_dir.
class RunAction final : public G4UserRunAction {
public:
    explicit RunAction(std::filesystem::path run_dir) : run_dir_(std::move(run_dir)) {}

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    [[nodiscard]] TruthRecorder& recorder() {
        return recorder_;
    }
    [[nodiscard]] DepositRecorder& deposits() {
        return deposits_;
    }
    [[nodiscard]] TruthTally& tally() {
        return tally_;
    }
    [[nodiscard]] NtupleWriter& writer() {
        return writer_;
    }

private:
    std::filesystem::path run_dir_;
    NtupleWriter writer_; // books the ntuple on construction, on the worker thread
    TruthRecorder recorder_;
    DepositRecorder deposits_;
    TruthTally tally_;
};

} // namespace riptide::sim
