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

#include <array>
#include <cstdint>
#include <string_view>

namespace riptide::schema {

inline constexpr std::int32_t kSchemaVersion = 1;
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

} // namespace riptide::schema
