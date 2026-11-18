// Compares the cross sections Geant4's high-precision neutron models read
// from G4NDL with the JEFF-3.3 evaluations G4NDL is built from. The checks
// of riptide_sim_check compute their expectations from JEFF-3.3 assuming
// that (i) the HP inelastic channel of carbon excludes radiative capture,
// which is a separate process, and (ii) Geant4 composes carbon from the
// natural-carbon file for C-12 plus the C-13 file. Both are pinned here, so
// a change of G4NDL that breaks them fails a test instead of biasing the
// checks.
//
// Without G4NDL (G4NEUTRONHPDATA unset) or the data disk every case is skipped.

#include <G4EnvironmentUtils.hh>
#include <G4ParticleHPManager.hh>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/endf_format.hpp"
#include "riptide/sim/endf_reader.hpp"
#include "riptide/sim/hp_data.hpp"

namespace {

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
namespace fs = std::filesystem;
namespace checks = riptide::sim::checks;
namespace endf = riptide::sim::endf;
namespace units = riptide::units;

constexpr const char* kHpDataEnv = "G4NEUTRONHPDATA";
constexpr const char* kHpCompressedSuffix = ".z";

constexpr int kCrossSectionMf = 3;
constexpr int kTotal = 1;
constexpr int kElastic = 2;
constexpr int kNonelastic = 3;
constexpr int kInelasticLevelsSum = 4;
constexpr int kCapture = 102;
// MT 600-849 are the level-resolved (n,p) ... (n,alpha) channels; when they
// are present MT 103-107 are their sums and must not be added again.
constexpr int kFirstChargedLevelMt = 600;
constexpr int kLastChargedLevelMt = 849;

constexpr double kEnMin = 10.0 * units::MeV;
constexpr double kEnMax = 20.0 * units::MeV;
const std::vector<double> kProbeEnergies = {10.0 * units::MeV, 14.0 * units::MeV,
                                            20.0 * units::MeV};
// Below the first inelastic threshold of C-12 (4.4 MeV level), where only
// elastic scattering and capture are open.
constexpr double kBelowThreshold = 1.0 * units::MeV;

// Elastic and capture of G4NDL are the JEFF-3.3 tables: observed equal to 0
// (H-1, C-0 elastic, C-13) and to 2.5e-6 (C-0 capture, rounding of the
// 7-digit values over a different grid).
constexpr double kCopiedRel = 1.0e-5;
// C-0 inelastic against MT3 - MT102: within 1e-7 at the probe energies,
// within 6.2e-4 at a few nodes near 10.5 MeV; the source of the latter is
// not known, the G4NDL nodes being a subset of the JEFF ones.
constexpr double kInelasticProbeRel = 1.0e-6;
constexpr double kInelasticNodeRel = 1.0e-3;
// Integral over 10-20 MeV: observed +2.5e-7 against MT3 - MT102 and
// -2.67e-4 against MT3, which is -sigma_capture / MT3. The two bounds sit
// far on either side of both values.
constexpr double kExcludesCaptureRel = 1.0e-5;
constexpr double kIncludesCaptureRel = 1.0e-4;
// C-13 inelastic is the sum of the JEFF-3.3 partial channels (observed
// 2.7e-4 at most), not MT3 - MT102 (off by up to 1.8e-2 near 14 MeV, since
// MT4 differs from the sum of MT51-91 in this evaluation).
constexpr double kC13PartialsRel = 5.0e-4;
constexpr double kC13SummaryOffMin = 1.0e-2;

std::optional<fs::path> hpDataDir() {
    const char* dir = G4FindDataDir(kHpDataEnv);
    if (dir == nullptr) {
        return std::nullopt;
    }
    return fs::path(dir);
}

fs::path hpDirOrSkip() {
    const std::optional<fs::path> dir = hpDataDir();
    if (!dir) {
        SKIP(std::string(kHpDataEnv) + " is not set: G4NDL not available");
    }
    // SKIP throws, so the optional is engaged here.
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return *dir;
}

// Channel ("Elastic", "Inelastic", "Capture") and nuclide file ("6_nat_Carbon").
fs::path hpFile(const fs::path& dir, const std::string& channel, const std::string& nuclide) {
    return dir / channel / "CrossSection" / (nuclide + kHpCompressedSuffix);
}

checks::HpCrossSection readHp(const std::string& channel, const std::string& nuclide) {
    const fs::path dir = hpDirOrSkip();
    const fs::path file = hpFile(dir, channel, nuclide);
    // GetDataStream would fall back to an uncompressed file or leave the
    // stream empty: require the file the HP model reads.
    REQUIRE(fs::is_regular_file(file));
    const std::string base = (dir / channel / "CrossSection" / nuclide).string();
    std::istringstream in;
    G4ParticleHPManager::GetInstance()->GetDataStream(base, in);
    return checks::parseHpCrossSection(in);
}

endf::EndfMaterial readJeffOrSkip(const std::string& name) {
    const fs::path path =
        riptide::data::resolveDataRoot(std::nullopt) / "reference" / "jeff33" / name;
    if (!fs::is_regular_file(path)) {
        SKIP("reference file not available: jeff33/" + name);
    }
    return endf::readEndfFile(path);
}

using Curve = std::function<double(double)>;

Curve jeffCurve(const endf::EndfMaterial& material, int mt) {
    endf::Tab1 table = endf::crossSection(material, mt);
    return [table = std::move(table)](double En) { return endf::evaluate(table, En); };
}

Curve hpCurve(const checks::HpCrossSection& table) {
    return [&table](double En) { return checks::evaluate(table, En); };
}

// Sorted, de-duplicated energies of the given grids within [kEnMin, kEnMax],
// ends included.
std::vector<double> unionGrid(const std::vector<std::vector<double>>& grids) {
    std::set<double> nodes = {kEnMin, kEnMax};
    for (const std::vector<double>& grid : grids) {
        for (const double En : grid) {
            if (En > kEnMin && En < kEnMax) {
                nodes.insert(En);
            }
        }
    }
    return {nodes.begin(), nodes.end()};
}

// Trapezoidal integral on the union grid: exact for lin-lin tables whose
// nodes are all in the grid.
double integral(const std::vector<double>& grid, const Curve& f) {
    double sum = 0.0;
    for (std::size_t i = 1; i < grid.size(); ++i) {
        sum += 0.5 * (f(grid[i - 1]) + f(grid[i])) * (grid[i] - grid[i - 1]);
    }
    return sum;
}

// HP nodes strictly inside [kEnMin, kEnMax].
std::vector<double> hpNodesInRange(const checks::HpCrossSection& table) {
    std::vector<double> nodes;
    std::copy_if(table.En.begin(), table.En.end(), std::back_inserter(nodes),
                 [](double En) { return En >= kEnMin && En <= kEnMax; });
    return nodes;
}

void checkCopied(const checks::HpCrossSection& hp, const endf::EndfMaterial& jeff, int mt) {
    const Curve reference = jeffCurve(jeff, mt);
    const std::vector<double> nodes = hpNodesInRange(hp);
    REQUIRE_FALSE(nodes.empty());
    for (const double En : nodes) {
        CAPTURE(mt, En);
        CHECK_THAT(checks::evaluate(hp, En), WithinRel(reference(En), kCopiedRel));
    }
}

bool isChargedLevelMt(int mt) {
    return mt >= kFirstChargedLevelMt && mt <= kLastChargedLevelMt;
}

// Every partial non-elastic channel of MF3: all reactions except the total,
// elastic, the sums MT3 and MT4, and capture.
std::vector<int> partialInelasticMts(const endf::EndfMaterial& material) {
    std::vector<int> mts;
    for (const auto& entry : material.sections) {
        const auto [mf, mt] = entry.first;
        if (mf != kCrossSectionMf) {
            continue;
        }
        REQUIRE_FALSE(isChargedLevelMt(mt));
        if (mt != kTotal && mt != kElastic && mt != kNonelastic && mt != kInelasticLevelsSum &&
            mt != kCapture) {
            mts.push_back(mt);
        }
    }
    return mts;
}

// Sum of the partials, each zero below its threshold (where its table starts).
Curve partialSum(const endf::EndfMaterial& material) {
    std::vector<endf::Tab1> tables;
    for (const int mt : partialInelasticMts(material)) {
        tables.push_back(endf::crossSection(material, mt));
    }
    return [tables = std::move(tables)](double En) {
        double sum = 0.0;
        for (const endf::Tab1& table : tables) {
            if (En >= table.x.front()) {
                sum += endf::evaluate(table, En);
            }
        }
        return sum;
    };
}

} // namespace

TEST_CASE("G4NDL has no C-12 file: the HP model uses natural carbon for it",
          "[hp_reference][carbon]") {
    const fs::path dir = hpDirOrSkip();
    for (const std::string channel : {"Elastic", "Inelastic", "Capture"}) {
        CAPTURE(channel);
        CHECK_FALSE(fs::exists(hpFile(dir, channel, "6_12_Carbon")));
        CHECK(fs::is_regular_file(hpFile(dir, channel, "6_nat_Carbon")));
        CHECK(fs::is_regular_file(hpFile(dir, channel, "6_13_Carbon")));
    }
}

TEST_CASE("G4NDL H-1 elastic and capture are the JEFF-3.3 tables", "[hp_reference][hydrogen]") {
    const endf::EndfMaterial jeff = readJeffOrSkip("n_0125_1-H-1.dat");
    checkCopied(readHp("Elastic", "1_1_Hydrogen"), jeff, kElastic);
    checkCopied(readHp("Capture", "1_1_Hydrogen"), jeff, kCapture);
}

TEST_CASE("G4NDL natural carbon elastic and capture are JEFF-3.3 C-0", "[hp_reference][carbon]") {
    const endf::EndfMaterial jeff = readJeffOrSkip("n_0600_6-C-0.dat");
    checkCopied(readHp("Elastic", "6_nat_Carbon"), jeff, kElastic);
    checkCopied(readHp("Capture", "6_nat_Carbon"), jeff, kCapture);
}

TEST_CASE("G4NDL natural carbon inelastic is JEFF-3.3 C-0 non-elastic without capture",
          "[hp_reference][carbon]") {
    const endf::EndfMaterial jeff = readJeffOrSkip("n_0600_6-C-0.dat");
    const checks::HpCrossSection hp = readHp("Inelastic", "6_nat_Carbon");
    const checks::HpCrossSection hp_capture = readHp("Capture", "6_nat_Carbon");
    const Curve nonelastic = jeffCurve(jeff, kNonelastic);
    const Curve capture = jeffCurve(jeff, kCapture);
    const Curve expected = [&](double En) { return nonelastic(En) - capture(En); };

    for (const double En : kProbeEnergies) {
        CAPTURE(En);
        CHECK_THAT(checks::evaluate(hp, En), WithinRel(expected(En), kInelasticProbeRel));
    }
    for (const double En : hpNodesInRange(hp)) {
        CAPTURE(En);
        CHECK_THAT(checks::evaluate(hp, En), WithinRel(expected(En), kInelasticNodeRel));
    }

    const std::vector<double> grid = unionGrid(
        {hp.En, endf::crossSection(jeff, kNonelastic).x, endf::crossSection(jeff, kCapture).x});
    const double hp_integral = integral(grid, hpCurve(hp));
    const double without_capture = integral(grid, expected);
    const double with_capture = integral(grid, nonelastic);
    CHECK(std::abs(hp_integral / without_capture - 1.0) < kExcludesCaptureRel);
    CHECK(std::abs(hp_integral / with_capture - 1.0) > kIncludesCaptureRel);

    // Below threshold the inelastic table has not started, capture has.
    CHECK(hp.En.front() > kBelowThreshold);
    CHECK_THAT(hp.sigma.front(), WithinAbs(0.0, 0.0));
    CHECK(checks::evaluate(hp_capture, kBelowThreshold) > 0.0);
}

TEST_CASE("G4NDL C-13 is JEFF-3.3 C-13 with inelastic from the partial channels",
          "[hp_reference][carbon]") {
    const endf::EndfMaterial jeff = readJeffOrSkip("n_0628_6-C-13.dat");
    checkCopied(readHp("Elastic", "6_13_Carbon"), jeff, kElastic);
    checkCopied(readHp("Capture", "6_13_Carbon"), jeff, kCapture);

    const checks::HpCrossSection hp = readHp("Inelastic", "6_13_Carbon");
    const Curve partials = partialSum(jeff);
    const Curve nonelastic = jeffCurve(jeff, kNonelastic);
    const Curve capture = jeffCurve(jeff, kCapture);
    double max_summary_off = 0.0;
    for (const double En : hpNodesInRange(hp)) {
        CAPTURE(En);
        const double value = checks::evaluate(hp, En);
        CHECK_THAT(value, WithinRel(partials(En), kC13PartialsRel));
        max_summary_off =
            std::max(max_summary_off, std::abs(value / (nonelastic(En) - capture(En)) - 1.0));
    }
    // Documents why the checks do not use MT3 - MT102 for C-13.
    CHECK(max_summary_off > kC13SummaryOffMin);
}
