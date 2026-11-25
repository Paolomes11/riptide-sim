#include "riptide/apps/case_b_geant4.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "riptide/io/run_info.hpp"

namespace riptide::apps {

namespace {

constexpr std::string_view kDirectionOption = "--direction";
constexpr char kSeparator = ',';
constexpr std::size_t kComponents = 3;

double norm(const geometry::Direction& v) noexcept {
    return std::sqrt(geometry::dot(v, v));
}

bool isUnit(const geometry::Direction& v, double tolerance) noexcept {
    return std::abs(norm(v) - 1.0) <= tolerance;
}

bool isFinite(const geometry::Direction& v) noexcept {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

std::vector<geometry::Direction> tracks(const std::vector<double>& x, const std::vector<double>& y,
                                        const std::vector<double>& z, const char* name) {
    std::vector<geometry::Direction> out;
    out.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        const geometry::Direction track{x[i], y[i], z[i]};
        if (!isUnit(track, kTrackNormTolerance)) {
            throw std::invalid_argument(std::string{"caseBSampleFromEvents: "} + name +
                                        " of event " + std::to_string(i) + " is not a unit vector");
        }
        out.push_back(track);
    }
    return out;
}

} // namespace

stats::CaseBPseudoSample caseBSampleFromEvents(const EventSample& events) {
    const std::size_t n = events.size();
    for (const auto* column :
         {&events.truth_En, &events.truth_Ep_phys, &events.truth_theta_p, &events.truth_p_dir_x,
          &events.truth_p_dir_y, &events.truth_p_dir_z, &events.reco_Ep, &events.reco_p_dir_x,
          &events.reco_p_dir_y, &events.reco_p_dir_z}) {
        if (column->size() != n) {
            throw std::invalid_argument(
                "caseBSampleFromEvents: a column does not have one entry per event");
        }
    }
    stats::CaseBPseudoSample sample;
    sample.truth_En = events.truth_En;
    sample.truth_Ep = events.truth_Ep_phys;
    sample.truth_theta_p = events.truth_theta_p;
    sample.truth_track =
        tracks(events.truth_p_dir_x, events.truth_p_dir_y, events.truth_p_dir_z, "truth_p_dir");
    sample.reco_Ep = events.reco_Ep;
    sample.reco_track =
        tracks(events.reco_p_dir_x, events.reco_p_dir_y, events.reco_p_dir_z, "reco_p_dir");
    return sample;
}

geometry::Direction parseDirection(std::string_view spec) {
    std::array<double, kComponents> values{};
    std::size_t start = 0;
    for (std::size_t c = 0; c < kComponents; ++c) {
        const std::size_t end = spec.find(kSeparator, start);
        const bool last = c + 1 == kComponents;
        if (last != (end == std::string_view::npos)) {
            throw CliError(std::string{kDirectionOption} + ": expected x,y,z, got '" +
                           std::string{spec} + "'");
        }
        const std::string field{spec.substr(start, last ? std::string_view::npos : end - start)};
        values.at(c) = parseCliNumber<double>(field, kDirectionOption);
        start = end + 1;
    }
    const geometry::Direction direction{values[0], values[1], values[2]};
    if (!isFinite(direction) || !isUnit(direction, kDirectionNormTolerance)) {
        throw CliError(std::string{kDirectionOption} + ": '" + std::string{spec} +
                       "' is not a unit vector");
    }
    return geometry::normalized(direction);
}

geometry::Direction fixedDirectionFromRunRecord(const nlohmann::json& record) {
    const std::string where = std::string{io::kRunRecordName} + ": config.source.direction";
    const nlohmann::json* value = nullptr;
    if (record.contains("config") && record.at("config").contains("source") &&
        record.at("config").at("source").contains("direction")) {
        value = &record.at("config").at("source").at("direction");
    }
    if (value == nullptr) {
        throw io::IoError(where + " is missing");
    }
    if (!value->is_array() || value->size() != kComponents) {
        throw io::IoError(where + " is not a fixed direction of three numbers");
    }
    std::array<double, kComponents> values{};
    for (std::size_t c = 0; c < kComponents; ++c) {
        if (!value->at(c).is_number()) {
            throw io::IoError(where + " is not a fixed direction of three numbers");
        }
        values.at(c) = value->at(c).get<double>();
    }
    const geometry::Direction direction{values[0], values[1], values[2]};
    if (!isFinite(direction) || !isUnit(direction, kDirectionNormTolerance)) {
        throw io::IoError(where + " is not a unit vector");
    }
    return geometry::normalized(direction);
}

void requireFixedDirection(const EventSample& events, const geometry::Direction& omega_n) {
    const std::size_t n = events.size();
    if (events.truth_omega_n_x.size() != n || events.truth_omega_n_y.size() != n ||
        events.truth_omega_n_z.size() != n) {
        throw std::invalid_argument(
            "requireFixedDirection: a column does not have one entry per event");
    }
    for (std::size_t i = 0; i < n; ++i) {
        const geometry::Direction omega{events.truth_omega_n_x[i], events.truth_omega_n_y[i],
                                        events.truth_omega_n_z[i]};
        // Written so that a NaN component fails the check.
        if (!(geometry::angleBetween(omega, omega_n) <= kFixedDirectionTolerance)) {
            throw std::invalid_argument("requireFixedDirection: event " + std::to_string(i) +
                                        " has a truth_omega_n off the configured direction");
        }
    }
}

} // namespace riptide::apps
