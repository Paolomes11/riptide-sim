#pragma once

// On-disk event schema: the contract between the simulation (writer) and the
// analysis (reader). Every column has a fixed name, storage type and unit in the
// project system (MeV, mm, ns, rad). Prefixes: truth_ = Monte Carlo truth,
// det_ = detector-level quantities, reco_ = reconstructed quantities.
//
// Any change of a column's name, type or meaning bumps kSchemaVersion, and
// readers must refuse files with a version they do not know. Only Int32 and
// Double columns exist because the Geant4 ntuple writer has no 64-bit integer
// column type.
//
// v1: first written schema (meta + truth columns).
// v2: adds kV2Branches (deposits, recoil-proton track, event class). The v1
//     columns keep name, type and meaning, so a v1 file is read as the first
//     kV1Columns columns of a v2 row.
// v3: adds kV3Branches (nuclear recoils made by the recoil proton through
//     electromagnetic Coulomb scattering). As for v2, older files are read as
//     the leading columns of a v3 row.
// v4: adds kV4Branches (start point of the primary neutron), which a
//     far-field source draws per event (D27). Older files are read as the
//     leading columns of a v4 row.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace riptide::schema {

inline constexpr std::int32_t kSchemaVersion = 4;
// Oldest version a reader of this build still accepts.
inline constexpr std::int32_t kOldestReadableSchemaVersion = 1;
inline constexpr std::string_view kEventTreeName = "events";

enum class ColumnType : std::uint8_t { Int32, Double };

struct BranchSpec {
    std::string_view name;
    ColumnType type;
    std::string_view unit; // "1" for dimensionless quantities and counters
    std::string_view description;
};

// Process of the neutron's first interaction inside the scintillator. Stored as
// Int32; the numeric values are part of the schema and must never be reused.
enum class FirstInteraction : std::int32_t {
    Unknown = 0,
    NpElastic = 1,
    NCElastic = 2,
    NCInelastic = 3,
    Other = 4,
};

// Class of an event from all hadronic interactions of the primary neutron
// inside the scintillator. Stored as Int32 like FirstInteraction: the numeric
// values are part of the schema and must never be reused. Precedence, highest
// first: any interaction outside the three below gives Other; then n-p with
// n-C gives Mixed; then the classes of a single kind.
enum class EventClass : std::int32_t {
    NoInteraction = 0,   // the primary neutron crossed without interacting
    SingleNp = 1,        // exactly one n-p elastic scatter, nothing else
    MultipleNp = 2,      // two or more n-p elastic scatters, nothing else
    NCElasticOnly = 3,   // n-C elastic scatters only
    NCInelasticNoNp = 4, // at least one n-C inelastic, no n-p (n-C elastic allowed)
    Mixed = 5,           // at least one n-p and at least one n-C
    Other = 6,           // at least one interaction not n-p or n-C elastic/inelastic
};

// Event identification. Multithreaded Geant4 runs write one file per thread in
// a non-deterministic order, so each event carries its id and the two seeds that
// reproduce it on its own. The seeds fit in 31 bits (see the simulation's
// seeding), so signed 32-bit columns hold them exactly.
inline constexpr std::array<BranchSpec, 3> kMetaBranches{{
    {"event_id", ColumnType::Int32, "1", "event number within the run"},
    {"rng_seed_0", ColumnType::Int32, "1", "first per-event random seed"},
    {"rng_seed_1", ColumnType::Int32, "1", "second per-event random seed"},
}};

inline constexpr std::array<BranchSpec, 16> kTruthBranches{{
    {"truth_En", ColumnType::Double, "MeV", "incident neutron kinetic energy"},
    {"truth_omega_n_x", ColumnType::Double, "1", "incident neutron direction, x"},
    {"truth_omega_n_y", ColumnType::Double, "1", "incident neutron direction, y"},
    {"truth_omega_n_z", ColumnType::Double, "1", "incident neutron direction, z"},
    {"truth_r0_x", ColumnType::Double, "mm", "first interaction vertex, x"},
    {"truth_r0_y", ColumnType::Double, "mm", "first interaction vertex, y"},
    {"truth_r0_z", ColumnType::Double, "mm", "first interaction vertex, z"},
    {"truth_t0", ColumnType::Double, "ns", "first interaction time"},
    {"truth_first_interaction", ColumnType::Int32, "1", "FirstInteraction code"},
    {"truth_n_np_scatters", ColumnType::Int32, "1", "number of n-p elastic scatters"},
    {"truth_Ep", ColumnType::Double, "MeV", "recoil proton kinetic energy (first n-p)"},
    {"truth_p_dir_x", ColumnType::Double, "1", "recoil proton direction, x"},
    {"truth_p_dir_y", ColumnType::Double, "1", "recoil proton direction, y"},
    {"truth_p_dir_z", ColumnType::Double, "1", "recoil proton direction, z"},
    {"truth_theta_p", ColumnType::Double, "rad", "proton recoil angle to the neutron direction"},
    {"truth_valid_np", ColumnType::Int32, "1",
     "1 if the truth_Ep/p_dir/theta_p columns are filled"},
}};

inline constexpr std::size_t kV1Columns = kMetaBranches.size() + kTruthBranches.size();

// Columns added in v2. The truth_np_vertex_* and truth_p_* columns describe
// the recoil proton of the first n-p scatter: like truth_Ep they are filled
// only when truth_valid_np is 1 and are 0 otherwise. The truth_p_* track
// quantities are the part of the proton track inside the scintillator.
inline constexpr std::array<BranchSpec, 18> kV2Branches{{
    {"det_Edep", ColumnType::Double, "MeV", "energy deposited in the scintillator, all particles"},
    {"det_Evis", ColumnType::Double, "MeV", "visible energy (Birks-quenched), all particles"},
    {"truth_np_vertex_x", ColumnType::Double, "mm", "first n-p scatter vertex, x"},
    {"truth_np_vertex_y", ColumnType::Double, "mm", "first n-p scatter vertex, y"},
    {"truth_np_vertex_z", ColumnType::Double, "mm", "first n-p scatter vertex, z"},
    {"truth_p_end_x", ColumnType::Double, "mm", "recoil proton end point in the scintillator, x"},
    {"truth_p_end_y", ColumnType::Double, "mm", "recoil proton end point in the scintillator, y"},
    {"truth_p_end_z", ColumnType::Double, "mm", "recoil proton end point in the scintillator, z"},
    {"truth_p_length", ColumnType::Double, "mm", "recoil proton path length in the scintillator"},
    {"truth_p_exits", ColumnType::Int32, "1", "1 if the recoil proton left the scintillator"},
    {"truth_p_Edep", ColumnType::Double, "MeV", "energy deposited by the recoil proton track"},
    {"truth_p_Evis", ColumnType::Double, "MeV", "visible energy of the recoil proton track"},
    {"truth_p_n_hadronic", ColumnType::Int32, "1",
     "hadronic interactions of the recoil proton in the scintillator"},
    {"truth_n_nC_elastic", ColumnType::Int32, "1", "number of n-C elastic scatters"},
    {"truth_n_nC_inelastic", ColumnType::Int32, "1", "number of n-C inelastic interactions"},
    {"truth_n_other", ColumnType::Int32, "1",
     "number of other hadronic interactions of the primary neutron"},
    {"truth_n_secondary_n_hadronic", ColumnType::Int32, "1",
     "hadronic interactions of non-primary neutrons in the scintillator"},
    {"truth_event_class", ColumnType::Int32, "1", "EventClass code"},
}};

inline constexpr std::size_t kV2Columns = kV1Columns + kV2Branches.size();

// Columns added in v3. Single Coulomb scattering of the recoil proton is an
// electromagnetic process, so truth_p_n_hadronic does not count it, yet it can
// give a nucleus (or a hydrogen proton) a recoil above threshold, which then
// carries energy away from the proton track. Like the other truth_p_* columns,
// filled only when truth_valid_np is 1.
inline constexpr std::array<BranchSpec, 1> kV3Branches{{
    {"truth_p_n_nuclear_recoil", ColumnType::Int32, "1",
     "non-hadronic steps of the recoil proton that made a nuclear recoil"},
}};

inline constexpr std::size_t kV3Columns = kV2Columns + kV3Branches.size();

// Columns added in v4. The checks need the chord of each neutron line through
// the scintillator, which with a far-field source changes from event to event
// and follows from the start point and truth_omega_n. Always filled; in a file
// older than v4 they read as 0, and the start point comes from the config.
inline constexpr std::array<BranchSpec, 3> kV4Branches{{
    {"truth_n_start_x", ColumnType::Double, "mm", "primary neutron start point, x"},
    {"truth_n_start_y", ColumnType::Double, "mm", "primary neutron start point, y"},
    {"truth_n_start_z", ColumnType::Double, "mm", "primary neutron start point, z"},
}};

inline constexpr std::size_t kV4Columns = kV3Columns + kV4Branches.size();

} // namespace riptide::schema
