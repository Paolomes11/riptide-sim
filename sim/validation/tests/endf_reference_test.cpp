// Reads the evaluated files on the data disk and compares a few numbers with
// those extracted by a separate, minimal parser written for this purpose.
// The disk may be absent (CI, unplugged SSD): then every case is skipped.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "riptide/common/data_paths.hpp"
#include "riptide/common/units.hpp"
#include "riptide/sim/endf_format.hpp"
#include "riptide/sim/endf_reader.hpp"

using Catch::Matchers::WithinRel;

namespace endf = riptide::sim::endf;
namespace units = riptide::units;

namespace {

// Values copied from the files are exact decimals, so the reader must give
// the same double; interpolated ones agree up to rounding in the two codes.
constexpr double kGridTolerance = 1e-12;
constexpr double kInterpolatedTolerance = 1e-9;
constexpr int kTotal = 1;
constexpr int kElastic = 2;
constexpr int kCapture = 102;

std::optional<std::filesystem::path> referenceFile(const std::string& library,
                                                   const std::string& name) {
    std::filesystem::path path =
        riptide::data::resolveDataRoot(std::nullopt) / "reference" / library / name;
    if (!std::filesystem::is_regular_file(path)) {
        return std::nullopt;
    }
    return path;
}

endf::EndfMaterial readOrSkip(const std::string& library, const std::string& name) {
    const std::optional<std::filesystem::path> path = referenceFile(library, name);
    if (!path) {
        SKIP("reference file not available: " + library + "/" + name);
    }
    // SKIP throws, so the optional is engaged here.
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return endf::readEndfFile(*path);
}

struct CrossSectionPoint {
    int mt;
    double En;    // MeV
    double barns; // value in the file
    bool on_grid; // En is a tabulated point
};

void checkCrossSections(const endf::EndfMaterial& material,
                        const std::vector<CrossSectionPoint>& points) {
    for (const CrossSectionPoint& point : points) {
        CAPTURE(point.mt, point.En);
        const double value = endf::evaluate(endf::crossSection(material, point.mt), point.En);
        CHECK_THAT(value / units::barn,
                   WithinRel(point.barns, point.on_grid ? kGridTolerance : kInterpolatedTolerance));
    }
}

struct LegendrePoint {
    double En; // MeV
    std::array<double, 4> a;
};

void checkLegendre(const endf::EndfMaterial& material, const std::vector<LegendrePoint>& points) {
    const endf::LegendreTable table = endf::legendreAngular(material);
    // The n-p coefficients are given in the centre-of-mass frame.
    CHECK(table.frame == endf::Frame::CenterOfMass);
    for (const LegendrePoint& point : points) {
        CAPTURE(point.En);
        const std::vector<double> a = endf::coefficients(table, point.En);
        REQUIRE(a.size() >= point.a.size());
        for (std::size_t i = 0; i < point.a.size(); ++i) {
            CAPTURE(i + 1);
            CHECK_THAT(a[i], WithinRel(point.a.at(i), kInterpolatedTolerance));
        }
    }
}

} // namespace

TEST_CASE("JEFF-3.3 H-1", "[reference][hydrogen]") {
    const endf::EndfMaterial material = readOrSkip("jeff33", "n_0125_1-H-1.dat");
    CHECK_THAT(material.za, WithinRel(1001.0, kGridTolerance));
    checkCrossSections(material, {{kElastic, 10.0, 0.9343974, true},
                                  {kElastic, 15.0, 0.6431586, true},
                                  {kElastic, 19.5, 0.495463, true},
                                  {kElastic, 20.0, 0.4827462, true},
                                  {kTotal, 10.0, 0.934429, true},
                                  {kTotal, 20.0, 0.4827735, true},
                                  {kCapture, 10.0, 3.162545e-05, true},
                                  {kCapture, 20.0, 2.722354e-05, true}});
    checkLegendre(material, {{10.0, {-1.153207e-02, 4.305945e-05, -2.609465e-04, 6.509103e-05}},
                             {15.0, {-1.527697e-02, 2.600421e-03, -7.874096e-04, 3.289702e-04}},
                             {19.5, {-1.609407e-02, 6.928997e-03, -1.403618e-03, 8.432785e-04}},
                             {20.0, {-1.608557e-02, 7.512273e-03, -1.468189e-03, 9.169295e-04}}});
}

TEST_CASE("ENDF/B-VIII.0 H-1", "[reference][hydrogen]") {
    const endf::EndfMaterial material = readOrSkip("endfb80", "n-001_H_001.endf");
    CHECK_THAT(material.za, WithinRel(1001.0, kGridTolerance));
    checkCrossSections(material, {{kElastic, 10.0, 0.9368747, true},
                                  {kElastic, 15.0, 0.6432359, true},
                                  {kElastic, 19.5, 0.4945948, true},
                                  {kElastic, 20.0, 0.4818408, true},
                                  {kTotal, 10.0, 0.9369061, true},
                                  {kTotal, 20.0, 0.4818679, true},
                                  {kCapture, 10.0, 3.143754e-05, true},
                                  {kCapture, 20.0, 2.710792e-05, true}});
    checkLegendre(material, {{10.0, {-1.064356e-02, 1.473377e-04, -1.442281e-04, 2.423436e-05}},
                             {15.0, {-1.506708e-02, 2.657640e-03, -3.815239e-04, 1.280747e-04}},
                             {19.5, {-1.616003e-02, 7.177100e-03, -7.276719e-04, 3.439591e-04}},
                             {20.0, {-1.612968e-02, 7.806923e-03, -7.769820e-04, 3.765550e-04}}});
}

TEST_CASE("JEFF-3.3 natural carbon and C-13", "[reference][carbon]") {
    const endf::EndfMaterial natural = readOrSkip("jeff33", "n_0600_6-C-0.dat");
    CHECK_THAT(natural.za, WithinRel(6000.0, kGridTolerance));
    checkCrossSections(natural, {{kTotal, 10.0, 1.172996, true},
                                 {kTotal, 20.0, 1.507006, true},
                                 {kElastic, 10.0, 0.6567711, true},
                                 {kElastic, 20.0, 0.9807161, true},
                                 {kCapture, 10.0, 8.47797e-05, true},
                                 {kCapture, 20.0, 0.000217352, true}});
    const endf::EndfMaterial c13 = readOrSkip("jeff33", "n_0628_6-C-13.dat");
    CHECK_THAT(c13.za, WithinRel(6013.0, kGridTolerance));
    checkCrossSections(c13, {{kTotal, 10.0, 1.539138, true},
                             {kTotal, 20.0, 1.527526, true},
                             {kElastic, 10.0, 0.826488, true},
                             {kElastic, 20.0, 0.8835171, true}});
}

TEST_CASE("ENDF/B-VIII.0 C-12 and C-13", "[reference][carbon]") {
    const endf::EndfMaterial c12 = readOrSkip("endfb80", "n-006_C_012.endf");
    CHECK_THAT(c12.za, WithinRel(6012.0, kGridTolerance));
    checkCrossSections(c12, {{kTotal, 10.0, 1.173, true},
                             {kTotal, 20.0, 1.50699997, true},
                             {kElastic, 10.0, 0.6567711, true},
                             {kElastic, 20.0, 0.9807161, true},
                             {kCapture, 10.0, 8.8896196e-05, false},
                             {kCapture, 20.0, 0.00021127, true}});
    const endf::EndfMaterial c13 = readOrSkip("endfb80", "n-006_C_013.endf");
    CHECK_THAT(c13.za, WithinRel(6013.0, kGridTolerance));
    checkCrossSections(c13, {{kTotal, 10.0, 1.17516089, true},
                             {kTotal, 20.0, 1.51767416, true},
                             {kElastic, 10.0, 1.032367, true},
                             {kElastic, 20.0, 1.373846, true}});
}
