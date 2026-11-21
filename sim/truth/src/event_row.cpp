#include "riptide/sim/event_row.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>

namespace riptide::sim {

const schema::BranchSpec& columnSpec(std::size_t i) {
    constexpr std::size_t n_meta = schema::kMetaBranches.size();
    if (i < n_meta) {
        return schema::kMetaBranches.at(i);
    }
    if (i < schema::kV1Columns) {
        return schema::kTruthBranches.at(i - n_meta);
    }
    if (i < schema::kV2Columns) {
        return schema::kV2Branches.at(i - schema::kV1Columns);
    }
    return schema::kV3Branches.at(i - schema::kV2Columns);
}

std::size_t columnCount(std::int32_t schema_version) {
    switch (schema_version) {
    case 1:
        return schema::kV1Columns;
    case 2:
        return schema::kV2Columns;
    case 3:
        return schema::kV3Columns;
    default:
        throw std::out_of_range("schema version " + std::to_string(schema_version) +
                                " is not readable by this build");
    }
}

std::size_t columnIndex(std::string_view name) {
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        if (columnSpec(i).name == name) {
            return i;
        }
    }
    throw std::out_of_range("no column named '" + std::string(name) + "' in the event schema");
}

namespace {

std::int32_t flag(bool b) {
    return b ? 1 : 0;
}

} // namespace

EventRow eventRow(const EventMeta& meta, const TruthRecord& record, const DepositRecord& deposit) {
    // Proton-side values pass through only for a valid n-p (D7).
    const bool np = record.valid_np;
    const Vec3 vertex = np ? record.np_vertex : Vec3{};
    const ProtonTrack p = np ? deposit.proton : ProtonTrack{};
    const auto event_class = static_cast<std::int32_t>(classifyEvent(interactionCounts(record)));
    // Same order as kMetaBranches, kTruthBranches, kV2Branches, kV3Branches; the unit test
    // checks every alternative against the schema's column type.
    return {
        meta.event_id,
        meta.rng_seed_0,
        meta.rng_seed_1,
        record.En,
        record.omega_n.x,
        record.omega_n.y,
        record.omega_n.z,
        record.r0.x,
        record.r0.y,
        record.r0.z,
        record.t0,
        static_cast<std::int32_t>(record.first_interaction),
        record.n_np_scatters,
        record.Ep,
        record.p_dir.x,
        record.p_dir.y,
        record.p_dir.z,
        record.theta_p,
        flag(record.valid_np),
        deposit.Edep,
        deposit.Evis,
        vertex.x,
        vertex.y,
        vertex.z,
        p.end.x,
        p.end.y,
        p.end.z,
        p.length,
        flag(p.exits),
        p.Edep,
        p.Evis,
        p.n_hadronic,
        record.n_nC_elastic,
        record.n_nC_inelastic,
        record.n_other,
        record.n_secondary_n_hadronic,
        event_class,
        p.n_nuclear_recoil,
    };
}

namespace {

double doubleColumn(const EventRow& row, std::string_view name) {
    return std::get<double>(row.at(columnIndex(name)));
}

std::int32_t intColumn(const EventRow& row, std::string_view name) {
    return std::get<std::int32_t>(row.at(columnIndex(name)));
}

bool flagColumn(const EventRow& row, std::string_view name) {
    const std::int32_t v = intColumn(row, name);
    if (v != 0 && v != 1) {
        throw std::out_of_range(std::string(name) + ": expected 0 or 1, got " + std::to_string(v));
    }
    return v == 1;
}

schema::FirstInteraction firstInteractionFromCode(std::int32_t code) {
    // Listed explicitly rather than range-checked, so a code added to the enum
    // without being handled here is not silently accepted.
    switch (code) {
    case static_cast<std::int32_t>(schema::FirstInteraction::Unknown):
        return schema::FirstInteraction::Unknown;
    case static_cast<std::int32_t>(schema::FirstInteraction::NpElastic):
        return schema::FirstInteraction::NpElastic;
    case static_cast<std::int32_t>(schema::FirstInteraction::NCElastic):
        return schema::FirstInteraction::NCElastic;
    case static_cast<std::int32_t>(schema::FirstInteraction::NCInelastic):
        return schema::FirstInteraction::NCInelastic;
    case static_cast<std::int32_t>(schema::FirstInteraction::Other):
        return schema::FirstInteraction::Other;
    default:
        throw std::out_of_range("truth_first_interaction: unknown code " + std::to_string(code));
    }
}

} // namespace

TruthRecord truthRecord(const EventRow& row) {
    TruthRecord r;
    r.En = doubleColumn(row, "truth_En");
    r.omega_n = {doubleColumn(row, "truth_omega_n_x"), doubleColumn(row, "truth_omega_n_y"),
                 doubleColumn(row, "truth_omega_n_z")};
    r.r0 = {doubleColumn(row, "truth_r0_x"), doubleColumn(row, "truth_r0_y"),
            doubleColumn(row, "truth_r0_z")};
    r.t0 = doubleColumn(row, "truth_t0");
    r.first_interaction = firstInteractionFromCode(intColumn(row, "truth_first_interaction"));
    r.n_np_scatters = intColumn(row, "truth_n_np_scatters");
    r.Ep = doubleColumn(row, "truth_Ep");
    r.p_dir = {doubleColumn(row, "truth_p_dir_x"), doubleColumn(row, "truth_p_dir_y"),
               doubleColumn(row, "truth_p_dir_z")};
    r.theta_p = doubleColumn(row, "truth_theta_p");
    r.valid_np = flagColumn(row, "truth_valid_np");
    r.np_vertex = {doubleColumn(row, "truth_np_vertex_x"), doubleColumn(row, "truth_np_vertex_y"),
                   doubleColumn(row, "truth_np_vertex_z")};
    r.n_nC_elastic = intColumn(row, "truth_n_nC_elastic");
    r.n_nC_inelastic = intColumn(row, "truth_n_nC_inelastic");
    r.n_other = intColumn(row, "truth_n_other");
    r.n_secondary_n_hadronic = intColumn(row, "truth_n_secondary_n_hadronic");
    return r;
}

DepositRecord depositRecord(const EventRow& row) {
    DepositRecord d;
    d.Edep = doubleColumn(row, "det_Edep");
    d.Evis = doubleColumn(row, "det_Evis");
    ProtonTrack& p = d.proton;
    p.seen = flagColumn(row, "truth_valid_np");
    p.start = {doubleColumn(row, "truth_np_vertex_x"), doubleColumn(row, "truth_np_vertex_y"),
               doubleColumn(row, "truth_np_vertex_z")};
    p.end = {doubleColumn(row, "truth_p_end_x"), doubleColumn(row, "truth_p_end_y"),
             doubleColumn(row, "truth_p_end_z")};
    p.length = doubleColumn(row, "truth_p_length");
    p.exits = flagColumn(row, "truth_p_exits");
    p.Edep = doubleColumn(row, "truth_p_Edep");
    p.Evis = doubleColumn(row, "truth_p_Evis");
    p.n_hadronic = intColumn(row, "truth_p_n_hadronic");
    p.n_nuclear_recoil = intColumn(row, "truth_p_n_nuclear_recoil");
    return d;
}

} // namespace riptide::sim
