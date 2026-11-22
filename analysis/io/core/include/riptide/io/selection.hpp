#pragma once

// Event selection as an ordered chain of named cuts. Every cut keeps its own
// count in the cut flow, so the effect of each one is visible and the chain
// written in analysis.json says exactly which events the analysis used.
//
// The cuts read truth columns on purpose: the first analysis compares the
// kinematics with the toy model on clean single n-p events, and the effect of
// relaxing each cut is measured later as a systematic.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace riptide::io {

enum class Cut {
    SingleNp,  // truth_event_class == SingleNp
    HasProton, // truth_valid_np == 1
    Contained, // truth_p_exits == 0
    EvisMin,   // det_Evis > threshold
};

[[nodiscard]] std::string_view cutName(Cut cut) noexcept;

struct CutSpec {
    Cut cut = Cut::SingleNp;
    double threshold = 0.0; // MeV; used by EvisMin only
};

// The columns the cuts read, for one event.
struct CutInputs {
    std::int32_t truth_event_class = 0;
    std::int32_t truth_p_exits = 0;
    double det_Evis = 0.0; // MeV
    // A single n-p event can still lack its proton track (the first n-p
    // produced no tracked proton); its truth Ep is then 0 and Case A cannot
    // describe its measured Ep, so the has_proton cut removes it.
    std::int32_t truth_valid_np = 0;
};

class Selection {
public:
    // Parses "single_np,has_proton,contained,evis_min=0.5": comma-separated cut names,
    // applied in the given order. evis_min needs an explicit threshold in MeV,
    // finite and >= 0. Throws IoError on an empty chain, an unknown or
    // repeated cut, or a malformed threshold.
    [[nodiscard]] static Selection parse(std::string_view text);

    [[nodiscard]] const std::vector<CutSpec>& cuts() const noexcept {
        return cuts_;
    }

    // Index of the first cut the event fails, or nullopt if it passes all.
    [[nodiscard]] std::optional<std::size_t> firstFailingCut(const CutInputs& event) const;

    // Canonical text, re-parsable to the same chain.
    [[nodiscard]] std::string toString() const;

private:
    std::vector<CutSpec> cuts_;
};

// Sequential cut flow: passed(k) counts the events that pass cuts 0..k.
// Counts can only decrease along the chain, by construction.
class CutFlow {
public:
    explicit CutFlow(std::size_t n_cuts) : passed_(n_cuts, 0) {}

    void record(std::optional<std::size_t> first_failing);

    [[nodiscard]] std::uint64_t total() const noexcept {
        return total_;
    }
    [[nodiscard]] std::uint64_t passed(std::size_t cut_index) const {
        return passed_.at(cut_index);
    }
    [[nodiscard]] std::size_t cuts() const noexcept {
        return passed_.size();
    }
    // Events that pass the whole chain.
    [[nodiscard]] std::uint64_t selected() const noexcept {
        return passed_.empty() ? total_ : passed_.back();
    }

private:
    std::uint64_t total_ = 0;
    std::vector<std::uint64_t> passed_;
};

} // namespace riptide::io
