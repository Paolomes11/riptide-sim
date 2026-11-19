#include "riptide/sim/endf_reader.hpp"

#include <algorithm>
#include <fstream>
#include <string_view>
#include <type_traits>
#include <utility>

#include "riptide/common/units.hpp"

namespace riptide::sim::endf {

namespace {

// MF and MT numbers used here (ENDF-102, Appendix B).
constexpr int kMfGeneral = 1;
constexpr int kMtDescription = 451;
constexpr int kMfCrossSection = 3;
constexpr int kMfAngular = 4;

// MF4 flags (ENDF-102, section 4.2).
constexpr int kLttLegendre = 1;
constexpr int kLiNotAllIsotropic = 0;

std::string sectionName(const EndfMaterial& material, int mf, int mt) {
    return material.source + " MF" + std::to_string(mf) + "/MT" + std::to_string(mt);
}

const EndfSection& findSection(const EndfMaterial& material, int mf, int mt) {
    const auto found = material.sections.find({mf, mt});
    if (found == material.sections.end()) {
        throw EndfError(sectionName(material, mf, mt) + ": section not in the file");
    }
    return found->second;
}

// The six fields of a CONT-like record (CONT, HEAD, LIST, TAB1, TAB2).
struct Control {
    double c1 = 0.0;
    double c2 = 0.0;
    int l1 = 0;
    int l2 = 0;
    int n1 = 0;
    int n2 = 0;
};

// Sequential reader of the records of one section; every error names the
// section and the line, so a bad file can be inspected by hand.
class RecordCursor {
public:
    RecordCursor(const EndfSection& section, std::string name)
        : section_(&section), name_(std::move(name)) {}

    Control control() {
        const EndfLine& line = next();
        constexpr std::size_t kC1 = 0;
        constexpr std::size_t kC2 = 1;
        constexpr std::size_t kL1 = 2;
        constexpr std::size_t kL2 = 3;
        constexpr std::size_t kN1 = 4;
        constexpr std::size_t kN2 = 5;
        return guarded([&] {
            return Control{line.real(kC1),    line.real(kC2),    line.integer(kL1),
                           line.integer(kL2), line.integer(kN1), line.integer(kN2)};
        });
    }

    // `count` values packed six per line, as in LIST and TAB1 bodies.
    std::vector<double> reals(std::size_t count) {
        std::vector<double> values;
        values.reserve(count);
        while (values.size() < count) {
            const EndfLine& line = next();
            for (std::size_t i = 0; i < kFieldsPerLine && values.size() < count; ++i) {
                values.push_back(guarded([&] { return line.real(i); }));
            }
        }
        return values;
    }

    std::vector<int> integers(std::size_t count) {
        std::vector<int> values;
        values.reserve(count);
        while (values.size() < count) {
            const EndfLine& line = next();
            for (std::size_t i = 0; i < kFieldsPerLine && values.size() < count; ++i) {
                values.push_back(guarded([&] { return line.integer(i); }));
            }
        }
        return values;
    }

    // NR (NBT, INT) pairs of a TAB1 or TAB2 record.
    void interpolationRanges(int ranges, std::vector<int>& boundaries,
                             std::vector<InterpolationLaw>& laws) {
        const std::vector<int> pairs = integers(2 * count(ranges, "NR"));
        boundaries.clear();
        laws.clear();
        for (std::size_t i = 0; i < pairs.size(); i += 2) {
            boundaries.push_back(pairs[i]);
            laws.push_back(guarded([&] { return interpolationLaw(pairs[i + 1]); }));
        }
    }

    // A TAB1 record; returns its control fields and fills `table`.
    Control tab1(Tab1& table) {
        const Control head = control();
        interpolationRanges(head.n1, table.boundaries, table.laws);
        const std::vector<double> pairs = reals(2 * count(head.n2, "NP"));
        table.x.clear();
        table.y.clear();
        for (std::size_t i = 0; i < pairs.size(); i += 2) {
            table.x.push_back(pairs[i]);
            table.y.push_back(pairs[i + 1]);
        }
        validate(table, where());
        return head;
    }

    // A record length from the file; `minimum` is 1 for most counts, 0 where
    // the format allows an empty list.
    [[nodiscard]] std::size_t count(int value, std::string_view what, int minimum = 1) const {
        if (value < minimum) {
            fail(std::string{what} + " = " + std::to_string(value) + " is below " +
                 std::to_string(minimum));
        }
        return static_cast<std::size_t>(value);
    }

    // Data after the last expected record mean the section was misread.
    void expectEnd() const {
        if (next_ != section_->lines.size()) {
            fail("unexpected records after the end of the data");
        }
    }

    [[noreturn]] void fail(const std::string& what) const {
        throw EndfError(where() + ": " + what);
    }

    [[nodiscard]] std::string where() const {
        // next_ is the number of lines already consumed; the last one is the
        // line being decoded.
        const std::size_t line = section_->first_line + (next_ == 0 ? 0 : next_ - 1);
        return name_ + " line " + std::to_string(line);
    }

private:
    const EndfLine& next() {
        if (next_ >= section_->lines.size()) {
            fail("section ends before the data it announces");
        }
        return section_->lines[next_++];
    }

    // Adds the position to errors from the field parsers.
    template <typename F> std::invoke_result_t<F> guarded(F&& read) const {
        try {
            return std::forward<F>(read)();
        } catch (const EndfError& error) {
            fail(error.what());
        }
    }

    const EndfSection* section_;
    std::string name_;
    std::size_t next_ = 0;
};

// End-of-section, file, material and tape records carry MT = 0; the tape
// identification line carries MF = MT = 0.
bool isEndRecord(const EndfLine& line) {
    return line.mt == 0;
}

} // namespace

EndfMaterial readEndf(std::istream& input, const std::string& source) {
    EndfMaterial material;
    material.source = source;
    std::pair<int, int> current{0, 0};
    std::string text;
    for (std::size_t number = 1; std::getline(input, text); ++number) {
        EndfLine line;
        try {
            line = parseEndfLine(text);
        } catch (const EndfError& error) {
            throw EndfError(source + " line " + std::to_string(number) + ": " + error.what());
        }
        if (isEndRecord(line)) {
            current = {0, 0};
            continue;
        }
        if (material.mat == 0) {
            material.mat = line.mat;
        } else if (line.mat != material.mat) {
            throw EndfError(source + " line " + std::to_string(number) + ": MAT " +
                            std::to_string(line.mat) + " after MAT " +
                            std::to_string(material.mat) +
                            "; only single-material files are supported");
        }
        const std::pair<int, int> key{line.mf, line.mt};
        if (key != current) {
            const auto [section, inserted] = material.sections.try_emplace(key);
            if (!inserted) {
                throw EndfError(source + " line " + std::to_string(number) + ": MF" +
                                std::to_string(line.mf) + "/MT" + std::to_string(line.mt) +
                                " appears twice");
            }
            section->second.first_line = number;
            current = key;
        }
        material.sections.at(key).lines.push_back(std::move(line));
    }
    if (input.bad()) {
        throw EndfError(source + ": read error");
    }
    const EndfSection& description = findSection(material, kMfGeneral, kMtDescription);
    RecordCursor cursor(description, sectionName(material, kMfGeneral, kMtDescription));
    const Control head = cursor.control();
    material.za = head.c1;
    material.awr = head.c2;
    return material;
}

EndfMaterial readEndfFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw EndfError("cannot open ENDF file " + path.string());
    }
    return readEndf(input, path.string());
}

Tab1 crossSection(const EndfMaterial& material, int mt) {
    // [MAT,3,MT/ ZA, AWR, 0, 0, 0, 0]HEAD, then [MAT,3,MT/ QM, QI, 0, LR, NR, NP/ E / sigma]TAB1.
    RecordCursor cursor(findSection(material, kMfCrossSection, mt),
                        sectionName(material, kMfCrossSection, mt));
    cursor.control();
    Tab1 table;
    cursor.tab1(table);
    cursor.expectEnd();
    // Energies in eV and cross sections in barn in the file.
    for (double& energy : table.x) {
        energy *= units::eV;
    }
    for (double& sigma : table.y) {
        sigma *= units::barn;
    }
    return table;
}

LegendreTable legendreAngular(const EndfMaterial& material, int mt) {
    // ENDF-102 section 4.2.1:
    // [MAT,4,MT/ ZA, AWR, 0, LTT, 0, 0]HEAD
    // [MAT,4,MT/ 0.0, AWR, LI, LCT, 0, 0]CONT
    // [MAT,4,MT/ 0.0, 0.0, 0, 0, NR, NE/ E_int]TAB2
    // NE x [MAT,4,MT/ T, E_i, LT, 0, NL, 0/ a_l(E_i)]LIST
    RecordCursor cursor(findSection(material, kMfAngular, mt),
                        sectionName(material, kMfAngular, mt));
    const Control head = cursor.control();
    if (head.l2 != kLttLegendre) {
        cursor.fail("LTT = " + std::to_string(head.l2) +
                    ": only Legendre coefficients over the whole range (LTT = 1) are supported");
    }
    const Control flags = cursor.control();
    if (flags.l1 != kLiNotAllIsotropic) {
        cursor.fail("LI = " + std::to_string(flags.l1) + " is not supported");
    }
    LegendreTable table;
    if (flags.l2 == static_cast<int>(Frame::Lab)) {
        table.frame = Frame::Lab;
    } else if (flags.l2 == static_cast<int>(Frame::CenterOfMass)) {
        table.frame = Frame::CenterOfMass;
    } else {
        cursor.fail("LCT = " + std::to_string(flags.l2) + " is neither 1 (LAB) nor 2 (CM)");
    }
    const Control grid = cursor.control();
    cursor.interpolationRanges(grid.n1, table.boundaries, table.laws);
    // The format requires lin-lin interpolation of the coefficients, which
    // keeps the interpolated distribution positive and allows negative a_l.
    for (const InterpolationLaw law : table.laws) {
        if (law != InterpolationLaw::LinLin) {
            cursor.fail("Legendre coefficients must be interpolated lin-lin (INT = 2)");
        }
    }
    const std::size_t energies = cursor.count(grid.n2, "NE");
    for (std::size_t i = 0; i < energies; ++i) {
        const Control list = cursor.control();
        table.energies.push_back(list.c2 * units::eV);
        table.coefficients.push_back(cursor.reals(cursor.count(list.n1, "NL", 0)));
    }
    cursor.expectEnd();
    validateGrid(table.boundaries, table.laws, table.energies,
                 sectionName(material, kMfAngular, mt));
    return table;
}

std::vector<double> coefficients(const LegendreTable& table, double En) {
    const Bracket where = bracket(table.boundaries, table.laws, table.energies, En);
    const std::vector<double>& lower = table.coefficients.at(where.lower);
    if (where.lower == where.upper) {
        return lower;
    }
    const std::vector<double>& upper = table.coefficients.at(where.upper);
    const std::size_t degree = std::max(lower.size(), upper.size());
    const auto at = [](const std::vector<double>& values, std::size_t index) {
        return index < values.size() ? values[index] : 0.0;
    };
    std::vector<double> result(degree);
    for (std::size_t l = 0; l < degree; ++l) {
        result[l] = interpolate(where.law, table.energies[where.lower], at(lower, l),
                                table.energies[where.upper], at(upper, l), En);
    }
    return result;
}

} // namespace riptide::sim::endf
