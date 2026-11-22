#include "riptide/io/response_summary.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "riptide/common/np_kinematics.hpp"
#include "riptide/io/run_info.hpp"

namespace riptide::io {

namespace {

void requireFilled(const std::vector<double>& column, std::size_t n, const char* name) {
    if (column.size() != n) {
        throw IoError(std::string{"column "} + name + " has not one entry per event");
    }
}

// Mean and unbiased variance of `values`, two passes for accuracy.
struct Moments {
    double mean = 0.0;
    double variance = 0.0;
    double max_abs = 0.0;
};

Moments moments(const std::vector<double>& values) {
    const auto n = static_cast<double>(values.size());
    Moments m;
    for (const double v : values) {
        m.mean += v;
        m.max_abs = std::max(m.max_abs, std::abs(v));
    }
    m.mean /= n;
    double sum_squares = 0.0;
    for (const double v : values) {
        sum_squares += (v - m.mean) * (v - m.mean);
    }
    m.variance = sum_squares / (n - 1.0);
    return m;
}

// Summary of N(0, sigma) residuals; the pulls are defined only for sigma > 0.
nlohmann::json gaussianResidualJson(const std::vector<double>& residuals, double sigma) {
    const auto m = moments(residuals);
    const auto n = static_cast<double>(residuals.size());
    nlohmann::json doc;
    doc["mean"] = m.mean;
    doc["std"] = std::sqrt(m.variance);
    doc["max_abs"] = m.max_abs;
    if (sigma > 0.0) {
        doc["pull_mean"] = m.mean / (sigma / std::sqrt(n));
        doc["pull_variance"] = (m.variance / (sigma * sigma) - 1.0) / std::sqrt(2.0 / (n - 1.0));
    } else {
        doc["pull_mean"] = nullptr;
        doc["pull_variance"] = nullptr;
    }
    return doc;
}

double dot(const Direction& a, const Direction& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Angle between unit vectors, accurate also for small angles (acos is not).
double angleBetween(const Direction& a, const Direction& b) noexcept {
    const Direction c{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    return std::atan2(std::sqrt(dot(c, c)), dot(a, b));
}

} // namespace

nlohmann::json hpConversionSummary(const EventSample& sample) {
    const auto n = sample.size();
    if (n == 0) {
        throw IoError("no events to summarise");
    }
    requireFilled(sample.truth_En, n, "truth_En");
    requireFilled(sample.truth_Ep, n, "truth_Ep");
    requireFilled(sample.truth_Ep_phys, n, "truth_Ep_phys");
    requireFilled(sample.truth_theta_p, n, "truth_theta_p");

    double offset_min = std::numeric_limits<double>::infinity();
    double offset_max = -std::numeric_limits<double>::infinity();
    double offset_sum = 0.0;
    std::uint64_t with_offset = 0;
    double lab_angle_max = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double Ep_phys = sample.truth_Ep_phys[i];
        if (!(Ep_phys > 0.0)) {
            continue;
        }
        const double offset = (sample.truth_Ep[i] - Ep_phys) / Ep_phys;
        offset_min = std::min(offset_min, offset);
        offset_max = std::max(offset_max, offset);
        offset_sum += offset;
        ++with_offset;
        const double Ep_lab = kinematics::protonEnergy(sample.truth_En[i], sample.truth_theta_p[i]);
        lab_angle_max = std::max(lab_angle_max, std::abs(Ep_phys - Ep_lab) / Ep_phys);
    }
    nlohmann::json doc;
    doc["events"] = n;
    doc["events_without_offset"] = n - with_offset;
    if (with_offset > 0) {
        doc["relative_offset"] = {{"min", offset_min},
                                  {"max", offset_max},
                                  {"mean", offset_sum / static_cast<double>(with_offset)}};
        doc["lab_angle_relative_max"] = lab_angle_max;
    } else {
        doc["relative_offset"] = nullptr;
        doc["lab_angle_relative_max"] = nullptr;
    }
    return doc;
}

nlohmann::json resolutionSummary(const EventSample& sample, const Resolution& resolution) {
    requireValidResolution(resolution);
    const auto n = sample.size();
    if (n < 2) {
        throw IoError("at least two events are needed for a variance");
    }
    requireFilled(sample.truth_Ep_phys, n, "truth_Ep_phys");
    requireFilled(sample.truth_theta_p, n, "truth_theta_p");
    requireFilled(sample.truth_p_dir_x, n, "truth_p_dir_x");
    requireFilled(sample.truth_p_dir_y, n, "truth_p_dir_y");
    requireFilled(sample.truth_p_dir_z, n, "truth_p_dir_z");
    requireFilled(sample.reco_Ep, n, "reco_Ep");
    requireFilled(sample.reco_theta_p, n, "reco_theta_p");
    requireFilled(sample.reco_p_dir_x, n, "reco_p_dir_x");
    requireFilled(sample.reco_p_dir_y, n, "reco_p_dir_y");
    requireFilled(sample.reco_p_dir_z, n, "reco_p_dir_z");

    std::vector<double> Ep_residuals(n);
    std::vector<double> theta_residuals(n);
    double scaled_sum = 0.0; // sum of delta^2 / (2 sigma^2)
    double angle_max = 0.0;
    double norm_deviation_max = 0.0;
    const double sigma = resolution.sigma_theta;
    for (std::size_t i = 0; i < n; ++i) {
        Ep_residuals[i] = sample.reco_Ep[i] - sample.truth_Ep_phys[i];
        theta_residuals[i] = sample.reco_theta_p[i] - sample.truth_theta_p[i];
        const Direction truth{sample.truth_p_dir_x[i], sample.truth_p_dir_y[i],
                              sample.truth_p_dir_z[i]};
        const Direction reco{sample.reco_p_dir_x[i], sample.reco_p_dir_y[i],
                             sample.reco_p_dir_z[i]};
        const double angle = angleBetween(truth, reco);
        angle_max = std::max(angle_max, angle);
        norm_deviation_max =
            std::max(norm_deviation_max, std::abs(std::sqrt(dot(reco, reco)) - 1.0));
        if (sigma > 0.0) {
            scaled_sum += angle * angle / (2.0 * sigma * sigma);
        }
    }

    nlohmann::json direction;
    direction["angle_max"] = angle_max;
    direction["norm_deviation_max"] = norm_deviation_max;
    if (sigma > 0.0) {
        const double scaled_mean = scaled_sum / static_cast<double>(n);
        direction["scaled_angle_squared_mean"] = scaled_mean;
        direction["pull_mean"] = (scaled_mean - 1.0) * std::sqrt(static_cast<double>(n));
    } else {
        direction["scaled_angle_squared_mean"] = nullptr;
        direction["pull_mean"] = nullptr;
    }

    nlohmann::json doc;
    doc["events"] = n;
    doc["Ep_MeV"] = gaussianResidualJson(Ep_residuals, resolution.sigma_Ep);
    doc["theta_p_rad"] = gaussianResidualJson(theta_residuals, resolution.sigma_theta);
    doc["p_dir"] = direction;
    return doc;
}

} // namespace riptide::io
