#include "riptide/sim/event_row.hpp"

#include <stdexcept>
#include <string>

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

} // namespace riptide::sim
