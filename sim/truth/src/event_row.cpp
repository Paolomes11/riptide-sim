#include "riptide/sim/event_row.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>

namespace riptide::sim {

const schema::BranchSpec& columnSpec(std::size_t i) {
    constexpr std::size_t n_meta = schema::kMetaBranches.size();
    return i < n_meta ? schema::kMetaBranches.at(i) : schema::kTruthBranches.at(i - n_meta);
}

std::size_t columnIndex(std::string_view name) {
    for (std::size_t i = 0; i < kEventColumns; ++i) {
        if (columnSpec(i).name == name) {
            return i;
        }
    }
    throw std::out_of_range("no column named '" + std::string(name) + "' in the event schema");
}

EventRow eventRow(const EventMeta& meta, const TruthRecord& record) {
    // Same order as kMetaBranches followed by kTruthBranches; the unit test
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
        static_cast<std::int32_t>(record.valid_np ? 1 : 0),
    };
}

namespace {

double doubleColumn(const EventRow& row, std::string_view name) {
    return std::get<double>(row.at(columnIndex(name)));
}

std::int32_t intColumn(const EventRow& row, std::string_view name) {
    return std::get<std::int32_t>(row.at(columnIndex(name)));
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
    const std::int32_t valid = intColumn(row, "truth_valid_np");
    if (valid != 0 && valid != 1) {
        throw std::out_of_range("truth_valid_np: expected 0 or 1, got " + std::to_string(valid));
    }
    r.valid_np = valid == 1;
    return r;
}

} // namespace riptide::sim
