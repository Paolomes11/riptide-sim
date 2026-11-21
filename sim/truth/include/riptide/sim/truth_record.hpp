#pragma once

// Monte Carlo truth of one event, built from the hadronic steps of the primary
// neutron inside the scintillator. Geant4-free on purpose: the Geant4 side only
// translates a G4Step into a HadronicStep, and all the bookkeeping below is
// unit-tested without running a simulation.

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>

#include "riptide/common/event_schema.hpp"
#include "riptide/sim/run_config.hpp"

namespace riptide::sim {

inline constexpr int kHydrogenZ = 1;
inline constexpr int kCarbonZ = 6;

// Category of a hadronic process, from its Geant4 process subtype. Process
// names are not used: they differ between physics lists (e.g. nCapture vs
// nCaptureHP).
enum class HadronicKind { Elastic, Inelastic, Capture, Fission, Other };

struct RecoilProton {
    double Ep = 0.0; // kinetic energy [MeV]
    Vec3 dir;        // unit momentum direction
};

// What the truth needs from one hadronic step of the primary neutron.
struct HadronicStep {
    HadronicKind kind = HadronicKind::Other;
    int target_Z = 0;  // atomic number of the struck nucleus
    Vec3 position;     // interaction point [mm]
    double time = 0.0; // global time of the interaction [ns]
    // Recoil proton among the step's secondaries, if any.
    std::optional<RecoilProton> proton;
};

// Event identification, mirroring the meta columns of event_schema.hpp.
struct EventMeta {
    std::int32_t event_id = 0;
    std::int32_t rng_seed_0 = 0;
    std::int32_t rng_seed_1 = 0;
};

// Mirrors the truth_* columns of event_schema.hpp, plus one diagnostic counter.
// The event class is not stored: it follows from the counters through
// classifyEvent, so it cannot disagree with them.
struct TruthRecord {
    double En = 0.0;
    Vec3 omega_n;
    Vec3 r0;
    double t0 = 0.0;
    schema::FirstInteraction first_interaction = schema::FirstInteraction::Unknown;
    std::int32_t n_np_scatters = 0;
    double Ep = 0.0;
    Vec3 p_dir;
    double theta_p = 0.0;
    bool valid_np = false;
    // Vertex of the first n-p scatter [mm]; it differs from r0 when an n-C
    // interaction came first. Meaningful only with valid_np.
    Vec3 np_vertex;
    std::int32_t n_nC_elastic = 0;
    std::int32_t n_nC_inelastic = 0;
    // Hadronic interactions of the primary outside the three kinds above
    // (capture, inelastic on hydrogen, ...).
    std::int32_t n_other = 0;
    // Hadronic interactions in the scintillator of neutrons other than the
    // primary. Geant4's high-precision inelastic models kill the primary and
    // emit the outgoing neutron as a new track, so its later interactions
    // (an n-p after n-C inelastic, say) are invisible to the counters above;
    // this counter says how often that happens.
    std::int32_t n_secondary_n_hadronic = 0;
    // Elastic scatters on hydrogen whose step had no recoil proton secondary.
    // Expected to stay zero; not written to the ntuple.
    std::int32_t n_np_without_proton = 0;
    // Start point of the primary neutron [mm].
    Vec3 n_start;
};

[[nodiscard]] schema::FirstInteraction classifyFirstInteraction(HadronicKind kind, int target_Z);

// Hadronic interactions of the primary neutron in one event, by kind.
struct InteractionCounts {
    std::int32_t np = 0;
    std::int32_t nC_elastic = 0;
    std::int32_t nC_inelastic = 0;
    std::int32_t other = 0;
};

[[nodiscard]] InteractionCounts interactionCounts(const TruthRecord& record);

// Event class with the precedence documented on schema::EventClass. Throws
// std::invalid_argument on a negative count, which no recorder can produce.
[[nodiscard]] schema::EventClass classifyEvent(const InteractionCounts& counts);

// Angle between two non-zero vectors [rad], in [0, pi]. The cosine is clamped
// to [-1, 1] so rounding on nearly parallel vectors cannot make acos return NaN.
[[nodiscard]] double thetaP(const Vec3& omega_n, const Vec3& p_dir);

// Accumulates the truth of one event. The first hadronic step fixes the first
// interaction; the first elastic scatter on hydrogen, even if preceded by
// other interactions, gives the recoil proton (theta_p is measured from the
// beam direction omega_n, which is what the detector reconstructs against).
class TruthRecorder {
public:
    // Starts a new event, discarding the previous one.
    void beginEvent(const EventMeta& meta, double En, const Vec3& omega_n, const Vec3& n_start);
    void onHadronicStep(const HadronicStep& step);
    // A hadronic interaction in the scintillator of a neutron that is not the
    // primary; see TruthRecord::n_secondary_n_hadronic.
    void onSecondaryNeutronHadronic();

    [[nodiscard]] const EventMeta& meta() const {
        return meta_;
    }
    [[nodiscard]] const TruthRecord& record() const {
        return record_;
    }

private:
    EventMeta meta_;
    TruthRecord record_;
};

// Per-thread run summary: how many events fell in each first-interaction
// category, and the diagnostic counter summed over events.
class TruthTally {
public:
    // One slot per FirstInteraction code, Unknown (no interaction) included.
    static constexpr std::size_t kCategories =
        static_cast<std::size_t>(schema::FirstInteraction::Other) + 1;

    void add(const TruthRecord& record);

    [[nodiscard]] std::int64_t events() const {
        return events_;
    }
    [[nodiscard]] std::int64_t count(schema::FirstInteraction category) const {
        return counts_.at(static_cast<std::size_t>(category));
    }
    [[nodiscard]] std::int64_t npWithoutProton() const {
        return np_without_proton_;
    }

private:
    std::int64_t events_ = 0;
    std::array<std::int64_t, kCategories> counts_{};
    std::int64_t np_without_proton_ = 0;
};

} // namespace riptide::sim
