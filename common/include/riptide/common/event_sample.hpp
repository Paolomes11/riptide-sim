#pragma once

// Selected events as contiguous columns (struct of arrays), the form in which
// the analysis receives them. It is filled by analysis/io and read by
// analysis/stats, which therefore needs neither ROOT nor the reader: this
// header is plain C++ on purpose.
//
// Names follow the event schema: truth_ columns are copied from the files,
// truth_Ep_phys is truth_Ep with the HP target-mass convention removed, reco_
// columns carry the detector resolution applied in the analysis.

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <vector>

namespace riptide {

struct EventSample {
    std::vector<std::int32_t> event_id;

    std::vector<double> truth_En;        // MeV
    std::vector<double> truth_omega_n_x; // 1
    std::vector<double> truth_omega_n_y; // 1
    std::vector<double> truth_omega_n_z; // 1

    std::vector<double> truth_Ep;      // MeV, as written by Geant4 (HP convention)
    std::vector<double> truth_Ep_phys; // MeV, physical n-p kinematics
    std::vector<double> truth_theta_p; // rad
    std::vector<double> truth_p_dir_x; // 1
    std::vector<double> truth_p_dir_y; // 1
    std::vector<double> truth_p_dir_z; // 1

    std::vector<double> det_Evis; // MeV

    std::vector<double> reco_Ep;      // MeV
    std::vector<double> reco_theta_p; // rad
    std::vector<double> reco_p_dir_x; // 1
    std::vector<double> reco_p_dir_y; // 1
    std::vector<double> reco_p_dir_z; // 1

    [[nodiscard]] std::size_t size() const noexcept {
        return event_id.size();
    }

    // Calls visit(column) on every double column, so checks and bulk operations
    // cannot forget one when a column is added.
    template <typename F> void forEachDoubleColumn(F&& visit) {
        forEach(*this, visit);
    }
    template <typename F> void forEachDoubleColumn(F&& visit) const {
        forEach(*this, visit);
    }

    // True when every column has one entry per event: a partly filled sample
    // is a bug, so a column not computed yet is not allowed to stay empty.
    [[nodiscard]] bool consistent() const {
        bool same = true;
        forEachDoubleColumn(
            [&](const std::vector<double>& column) { same = same && column.size() == size(); });
        return same;
    }

    void reserve(std::size_t n) {
        event_id.reserve(n);
        forEachDoubleColumn([n](std::vector<double>& column) { column.reserve(n); });
    }

private:
    // One list of the columns for the const and non-const visitors.
    template <typename Self, typename F> static void forEach(Self& self, F& visit) {
        for (auto* column :
             {&self.truth_En, &self.truth_omega_n_x, &self.truth_omega_n_y, &self.truth_omega_n_z,
              &self.truth_Ep, &self.truth_Ep_phys, &self.truth_theta_p, &self.truth_p_dir_x,
              &self.truth_p_dir_y, &self.truth_p_dir_z, &self.det_Evis, &self.reco_Ep,
              &self.reco_theta_p, &self.reco_p_dir_x, &self.reco_p_dir_y, &self.reco_p_dir_z}) {
            visit(*column);
        }
    }
};

} // namespace riptide
