#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstddef>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/common/units.hpp"
#include "riptide/sim/endf_reader.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace endf = riptide::sim::endf;
namespace units = riptide::units;
using endf::EndfError;

namespace {

constexpr int kMat = 9001; // fictitious material number

// Right-justifies a value in an 11-column field.
std::string field(std::string_view value) {
    return std::string(endf::kFieldWidth - value.size(), ' ') + std::string{value};
}

std::string field(int value) {
    return field(std::to_string(value));
}

std::string tail(int mat, int mf, int mt) {
    std::string mat_text = std::to_string(mat);
    std::string mf_text = std::to_string(mf);
    std::string mt_text = std::to_string(mt);
    return std::string(4 - mat_text.size(), ' ') + mat_text + std::string(2 - mf_text.size(), ' ') +
           mf_text + std::string(3 - mt_text.size(), ' ') + mt_text;
}

// One record from up to six already formatted fields.
std::string record(const std::vector<std::string>& fields, int mat, int mf, int mt) {
    std::string data;
    for (const std::string& value : fields) {
        data += value;
    }
    data.resize(endf::kDataWidth, ' ');
    return data + tail(mat, mf, mt);
}

struct TapeOptions {
    int ltt = 1;               // MF4 representation flag
    int legendre_law = 2;      // INT of the MF4 TAB2 record
    int mf3_mat = kMat;        // MAT of the MF3 records
    bool repeat_mf3 = false;   // MF3/MT2 given a second time after MF4
    bool truncate_mf4 = false; // last LIST body line missing
    bool extra_record = false; // spurious record at the end of MF3/MT2
};

// A minimal single-material tape: description, MF3/MT2 with two
// interpolation ranges, MF4/MT2 with three energies and NL = 1, 3, 2.
std::string syntheticTape(const TapeOptions& options = {}) {
    const std::string zero = field("0.000000+0");
    const std::string za = field("1.001000+3");
    const std::string awr = field("9.991673-1");
    std::vector<std::string> lines;
    const auto end_section = [&](int mf) { lines.push_back(record({}, kMat, mf, 0)); };
    const auto end_file = [&] { lines.push_back(record({}, kMat, 0, 0)); };

    lines.push_back(record({" synthetic tape for the reader tests"}, 1, 0, 0));
    lines.push_back(record({za, awr, field(0), field(0), field(0), field(1)}, kMat, 1, 451));
    lines.push_back(record({" free text, not numbers: 1.0-2 abc"}, kMat, 1, 451));
    end_section(1);
    end_file();

    // sigma = 8 b * (E / 1 MeV)^(-1/2) log-log up to 16 MeV, then lin-lin.
    const auto mf3 = [&] {
        const int mat = options.mf3_mat;
        lines.push_back(record({za, awr, field(0), field(0), field(0), field(0)}, mat, 3, 2));
        lines.push_back(record({zero, zero, field(0), field(0), field(2), field(5)}, mat, 3, 2));
        lines.push_back(record({field(3), field(5), field(5), field(2)}, mat, 3, 2));
        lines.push_back(record({field("1.000000+6"), field("8.000000+0"), field("4.000000+6"),
                                field("4.000000+0"), field("1.600000+7"), field("2.000000+0")},
                               mat, 3, 2));
        lines.push_back(record(
            {field("2.000000+7"), field("6.000000+0"), field("3.000000+7"), field("1.000000+0")},
            mat, 3, 2));
        if (options.extra_record) {
            lines.push_back(record({zero, zero}, mat, 3, 2));
        }
        end_section(3);
    };
    mf3();
    end_file();

    lines.push_back(
        record({za, awr, field(0), field(options.ltt), field(0), field(0)}, kMat, 4, 2));
    lines.push_back(record({zero, awr, field(0), field(2), field(0), field(0)}, kMat, 4, 2));
    lines.push_back(record({zero, zero, field(0), field(0), field(1), field(3)}, kMat, 4, 2));
    lines.push_back(record({field(3), field(options.legendre_law)}, kMat, 4, 2));
    lines.push_back(
        record({zero, field("1.000000+6"), field(0), field(0), field(1), field(0)}, kMat, 4, 2));
    lines.push_back(record({field("1.000000-2")}, kMat, 4, 2));
    lines.push_back(
        record({zero, field("1.000000+7"), field(0), field(0), field(3), field(0)}, kMat, 4, 2));
    lines.push_back(
        record({field("2.000000-2"), field("-4.00000-3"), field("1.000000-3")}, kMat, 4, 2));
    lines.push_back(
        record({zero, field("2.000000+7"), field(0), field(0), field(2), field(0)}, kMat, 4, 2));
    if (!options.truncate_mf4) {
        lines.push_back(record({field("4.000000-2"), field("-6.00000-3")}, kMat, 4, 2));
    }
    end_section(4);
    end_file();

    if (options.repeat_mf3) {
        mf3();
        end_file();
    }
    lines.push_back(record({}, 0, 0, 0));  // MEND
    lines.push_back(record({}, -1, 0, 0)); // TEND

    std::string tape;
    for (const std::string& line : lines) {
        tape += line + '\n';
    }
    return tape;
}

template <typename Change> TapeOptions with(Change change) {
    TapeOptions options;
    change(options);
    return options;
}

endf::EndfMaterial read(const std::string& tape) {
    std::istringstream input(tape);
    return endf::readEndf(input, "synthetic");
}

} // namespace

TEST_CASE("a tape is grouped into sections of one material", "[endf][reader]") {
    const endf::EndfMaterial material = read(syntheticTape());
    CHECK(material.mat == kMat);
    CHECK_THAT(material.za, WithinULP(1001.0, 0));
    CHECK_THAT(material.awr, WithinULP(9.991673e-1, 0));
    REQUIRE(material.sections.size() == 3);
    CHECK(material.sections.at({3, 2}).lines.size() == 5);
    CHECK(material.sections.at({3, 2}).first_line == 6);
    CHECK(material.sections.at({4, 2}).lines.size() == 10);
}

TEST_CASE("MF3 cross section in MeV and mm^2", "[endf][reader]") {
    const endf::Tab1 sigma = endf::crossSection(read(syntheticTape()), 2);
    REQUIRE(sigma.x.size() == 5);
    CHECK_THAT(sigma.x.front(), WithinULP(1.0 * units::MeV, 0));
    // The upper end of the RIPTIDE domain must convert exactly, or a check at
    // 20 MeV would fall outside the table.
    CHECK_THAT(sigma.x[3], WithinULP(20.0 * units::MeV, 0));
    CHECK_THAT(endf::evaluate(sigma, 20.0 * units::MeV), WithinRel(6.0 * units::barn, 1e-15));
    // Log-log range: 8 b / sqrt(9).
    CHECK_THAT(endf::evaluate(sigma, 9.0 * units::MeV), WithinRel(8.0 / 3.0 * units::barn, 1e-13));
    // Lin-lin range.
    CHECK_THAT(endf::evaluate(sigma, 25.0 * units::MeV), WithinRel(3.5 * units::barn, 1e-14));
    CHECK_THROWS_WITH(endf::crossSection(read(syntheticTape()), 102),
                      ContainsSubstring("MF3/MT102") && ContainsSubstring("not in the file"));
}

TEST_CASE("MF4 Legendre coefficients", "[endf][reader]") {
    const endf::LegendreTable table = endf::legendreAngular(read(syntheticTape()));
    CHECK(table.frame == endf::Frame::CenterOfMass);
    REQUIRE(table.energies.size() == 3);
    CHECK_THAT(table.energies[1], WithinULP(10.0 * units::MeV, 0));
    CHECK(table.coefficients[0].size() == 1);
    CHECK(table.coefficients[1].size() == 3);
    CHECK(table.coefficients[2].size() == 2);

    SECTION("tabulated energy") {
        const std::vector<double> a = endf::coefficients(table, 20.0 * units::MeV);
        REQUIRE(a.size() == 2);
        CHECK_THAT(a[0], WithinULP(4.0e-2, 0));
        CHECK_THAT(a[1], WithinULP(-6.0e-3, 0));
    }
    SECTION("between NL = 1 and NL = 3: missing coefficients are zero") {
        const std::vector<double> a = endf::coefficients(table, 5.5 * units::MeV);
        REQUIRE(a.size() == 3);
        CHECK_THAT(a[0], WithinRel(1.5e-2, 1e-14));
        CHECK_THAT(a[1], WithinRel(-2.0e-3, 1e-14));
        CHECK_THAT(a[2], WithinRel(5.0e-4, 1e-14));
    }
    SECTION("between NL = 3 and NL = 2") {
        const std::vector<double> a = endf::coefficients(table, 15.0 * units::MeV);
        REQUIRE(a.size() == 3);
        CHECK_THAT(a[0], WithinRel(3.0e-2, 1e-14));
        CHECK_THAT(a[1], WithinRel(-5.0e-3, 1e-14));
        CHECK_THAT(a[2], WithinRel(5.0e-4, 1e-13));
    }
    SECTION("outside the energy grid") {
        CHECK_THROWS_AS(endf::coefficients(table, 0.5 * units::MeV), EndfError);
        CHECK_THROWS_AS(endf::coefficients(table, 20.5 * units::MeV), EndfError);
    }
}

TEST_CASE("unsupported or inconsistent files are refused", "[endf][reader]") {
    SECTION("LTT = 3 (Legendre plus tabulated)") {
        CHECK_THROWS_WITH(
            endf::legendreAngular(read(syntheticTape(with([](TapeOptions& o) { o.ltt = 3; })))),
            ContainsSubstring("LTT = 3"));
    }
    SECTION("Legendre coefficients not interpolated lin-lin") {
        CHECK_THROWS_WITH(endf::legendreAngular(read(
                              syntheticTape(with([](TapeOptions& o) { o.legendre_law = 5; })))),
                          ContainsSubstring("INT = 2"));
    }
    SECTION("second material") {
        CHECK_THROWS_WITH(read(syntheticTape(with([](TapeOptions& o) { o.mf3_mat = kMat + 1; }))),
                          ContainsSubstring("line 6") && ContainsSubstring("single-material"));
    }
    SECTION("section split in two") {
        CHECK_THROWS_WITH(read(syntheticTape(with([](TapeOptions& o) { o.repeat_mf3 = true; }))),
                          ContainsSubstring("MF3/MT2 appears twice"));
    }
    SECTION("truncated section") {
        CHECK_THROWS_WITH(endf::legendreAngular(read(
                              syntheticTape(with([](TapeOptions& o) { o.truncate_mf4 = true; })))),
                          ContainsSubstring("MF4/MT2") && ContainsSubstring("ends before"));
    }
    SECTION("records left over after the data") {
        CHECK_THROWS_WITH(
            endf::crossSection(
                read(syntheticTape(with([](TapeOptions& o) { o.extra_record = true; }))), 2),
            ContainsSubstring("unexpected records"));
    }
    SECTION("short line") {
        std::string tape = syntheticTape();
        const std::size_t third_line = tape.find('\n', tape.find('\n') + 1) + 1;
        tape.erase(third_line + 70, 5); // line 3 cut to 70 columns
        CHECK_THROWS_WITH(read(tape), ContainsSubstring("synthetic line 3"));
    }
    SECTION("no description section") {
        CHECK_THROWS_WITH(read(record({}, -1, 0, 0) + "\n"), ContainsSubstring("MF1/MT451"));
    }
}
