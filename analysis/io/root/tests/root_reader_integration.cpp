// Integration test of readEventColumns on the runs of the riptide_sim
// integration test (driven by run_io_root_it.sh). These check that the reader
// returns every event of a run once, in event_id order, independently of how
// the events were spread over files and threads, and that it refuses event
// files without a schema column.

#include <ROOT/RDataFrame.hxx>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/np_kinematics.hpp"
#include "riptide/io/event_columns.hpp"
#include "riptide/io/response.hpp"
#include "riptide/io/root_reader.hpp"
#include "riptide/io/run_info.hpp"
#include "riptide/io/selection.hpp"

namespace {

using riptide::io::EventColumns;
using riptide::io::IoError;
using riptide::io::readEventColumns;
using riptide::io::RunInfo;

constexpr std::string_view kRun1t = "ctest_it_1t";
constexpr std::string_view kRun2t = "ctest_it_2t";
constexpr std::uint64_t kBoundedEvents = 100;
// Relative HP offset (Ep_hp - Ep_phys) / Ep_phys expected over 10-19.5 MeV
// (phase-3 predictions), and the agreement of the physical energy with the
// bare-proton kinematics at the same lab angle; as in response_test.
constexpr double kOffsetMin = 5.33500e-4;
constexpr double kOffsetMax = 5.50179e-4;
constexpr double kLabAngleTolerance = 6e-6;

std::filesystem::path pathFromEnv(const char* name) {
    const char* value = std::getenv(name); // NOLINT(concurrency-mt-unsafe): single-threaded
    if (value == nullptr || *value == '\0') {
        throw std::runtime_error(std::string(name) + " is not set; run through run_io_root_it.sh");
    }
    return value;
}

RunInfo runInfo(std::string_view run_id) {
    return riptide::io::readRunInfo(pathFromEnv("RIPTIDE_IT_DATA_ROOT"), run_id, std::nullopt);
}

// Bit-for-bit equality of two double columns, without a float comparison.
bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

void checkIdentical(EventColumns& a, EventColumns& b) {
    CHECK(a.event_id == b.event_id);
    CHECK(a.truth_event_class == b.truth_event_class);
    CHECK(a.truth_p_exits == b.truth_p_exits);
    std::vector<const std::vector<double>*> a_doubles;
    a.forEachDoubleColumn(
        [&](std::string_view, std::vector<double>& column) { a_doubles.push_back(&column); });
    std::size_t index = 0;
    b.forEachDoubleColumn([&](std::string_view name, std::vector<double>& column) {
        INFO("column " << name);
        CHECK(sameBits(*a_doubles[index++], column));
    });
}

// Writes the first run's events to `file` through `frame_edit`, which drops or
// retypes a column, and returns a RunInfo that points at it.
template <typename Edit> RunInfo runWithEditedFile(const std::string& file, Edit frame_edit) {
    auto run = runInfo(kRun1t);
    std::vector<std::string> inputs;
    for (const auto& path : run.event_files) {
        inputs.push_back(path.string());
    }
    ROOT::RDataFrame frame(std::string{riptide::schema::kEventTreeName}, inputs);
    frame_edit(frame, file);
    run.event_files = {pathFromEnv("RIPTIDE_IT_SCRATCH") / file};
    return run;
}

} // namespace

TEST_CASE("every event is read once in event_id order", "[io_root]") {
    const auto run = runInfo(kRun1t);
    const auto columns = readEventColumns(run, 0);
    REQUIRE(columns.size() == run.events);
    for (std::size_t i = 0; i < columns.size(); ++i) {
        REQUIRE(columns.event_id[i] == static_cast<std::int32_t>(i));
    }
}

TEST_CASE("one and two threads give the same columns", "[io_root]") {
    const auto run_2t = runInfo(kRun2t);
    // The test is only meaningful if the events really are spread over files.
    REQUIRE(run_2t.event_files.size() == 2);
    auto columns_1t = readEventColumns(runInfo(kRun1t), 0);
    auto columns_2t = readEventColumns(run_2t, 0);
    checkIdentical(columns_1t, columns_2t);
}

TEST_CASE("max_events keeps the first events", "[io_root]") {
    auto all = readEventColumns(runInfo(kRun2t), 0);
    auto bounded = readEventColumns(runInfo(kRun2t), kBoundedEvents);
    REQUIRE(bounded.size() == kBoundedEvents);
    all.forEachInt32Column(
        [](std::string_view, std::vector<std::int32_t>& column) { column.resize(kBoundedEvents); });
    all.forEachDoubleColumn(
        [](std::string_view, std::vector<double>& column) { column.resize(kBoundedEvents); });
    checkIdentical(all, bounded);
}

TEST_CASE("the HP conversion of read events removes the expected offset", "[io_root]") {
    const auto run = runInfo(kRun1t);
    const auto columns = readEventColumns(run, 0);
    const auto selection = riptide::io::Selection::parse("single_np");
    riptide::io::CutFlow flow(selection.cuts().size());
    auto sample = riptide::io::selectEvents(columns, selection, flow);
    REQUIRE(sample.size() > 0);
    riptide::io::convertHpEnergies(sample, run.hp.target_mass_MeV);
    for (std::size_t i = 0; i < sample.size(); ++i) {
        INFO("event " << sample.event_id[i]);
        const double offset =
            (sample.truth_Ep[i] - sample.truth_Ep_phys[i]) / sample.truth_Ep_phys[i];
        CHECK(offset > kOffsetMin);
        CHECK(offset < kOffsetMax);
        CHECK_THAT(sample.truth_Ep_phys[i],
                   Catch::Matchers::WithinRel(riptide::kinematics::protonEnergy(
                                                  sample.truth_En[i], sample.truth_theta_p[i]),
                                              kLabAngleTolerance));
    }
}

TEST_CASE("a missing or retyped column is refused", "[io_root]") {
    const auto scratch = pathFromEnv("RIPTIDE_IT_SCRATCH");
    const auto missing = runWithEditedFile(
        "missing_det_Evis.root", [&](ROOT::RDataFrame& frame, const std::string& file) {
            std::vector<std::string> kept;
            for (const auto& name : frame.GetColumnNames()) {
                if (name != "det_Evis") {
                    kept.push_back(name);
                }
            }
            frame.Snapshot(std::string{riptide::schema::kEventTreeName}, (scratch / file).string(),
                           kept);
        });
    CHECK_THROWS_AS(readEventColumns(missing, 0), IoError);

    const auto retyped = runWithEditedFile(
        "float_det_Evis.root", [&](ROOT::RDataFrame& frame, const std::string& file) {
            frame
                .Redefine("det_Evis", [](double value) { return static_cast<float>(value); },
                          {"det_Evis"})
                .Snapshot(std::string{riptide::schema::kEventTreeName}, (scratch / file).string());
        });
    CHECK_THROWS_AS(readEventColumns(retyped, 0), IoError);
}
