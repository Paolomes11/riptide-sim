#include <catch2/catch_test_macros.hpp>

#include <stdlib.h> // mkdtemp (POSIX)

#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>

#include "riptide/common/data_paths.hpp"

namespace fs = std::filesystem;
using riptide::data::DataPathError;

namespace {

// Fresh directory under the system temp dir, removed at scope exit.
class TempDir {
public:
    TempDir() {
        // mkdtemp needs a mutable, NUL-terminated template.
        std::string tmpl = (fs::temp_directory_path() / "riptide-test-XXXXXX").string();
        if (::mkdtemp(tmpl.data()) == nullptr) {
            throw std::runtime_error("mkdtemp failed");
        }
        path_ = tmpl;
    }
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(path_, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&&) = delete;
    TempDir& operator=(TempDir&&) = delete;
    [[nodiscard]] const fs::path& path() const {
        return path_;
    }

private:
    fs::path path_;
};

} // namespace

TEST_CASE("data root precedence: cli > env > default", "[data_paths]") {
    using riptide::data::resolveDataRoot;
    const std::optional<std::string> none;
    CHECK(resolveDataRoot(std::string{"/a/cli"}, "/b/env") == fs::path{"/a/cli"});
    CHECK(resolveDataRoot(none, "/b/env") == fs::path{"/b/env"});
    CHECK(resolveDataRoot(none, nullptr) == fs::path{std::string{riptide::data::kDefaultDataRoot}});
    // Empty values count as unset.
    CHECK(resolveDataRoot(std::string{}, "") ==
          fs::path{std::string{riptide::data::kDefaultDataRoot}});
    CHECK(resolveDataRoot(none, "/b/env/../env2/") == fs::path{"/b/env2/"});
}

TEST_CASE("relative data roots are rejected", "[data_paths]") {
    using riptide::data::resolveDataRoot;
    CHECK_THROWS_AS(resolveDataRoot(std::string{"data"}, nullptr), DataPathError);
    CHECK_THROWS_AS(resolveDataRoot(std::nullopt, "./data"), DataPathError);
}

TEST_CASE("a missing root is an error and is not created", "[data_paths]") {
    const TempDir tmp;
    const fs::path missing = tmp.path() / "not-mounted";
    CHECK_THROWS_AS(riptide::data::requireWritableRoot(missing, 0, 0, true), DataPathError);
    CHECK_FALSE(fs::exists(missing));
}

TEST_CASE("the system disk is refused unless allowed", "[data_paths]") {
    const TempDir tmp;
    if (!riptide::data::sameDevice(tmp.path(), "/")) {
        SKIP("temp dir is not on the system disk on this machine");
    }
    CHECK_THROWS_AS(riptide::data::requireWritableRoot(tmp.path(), 0, 0, false), DataPathError);
    CHECK_NOTHROW(riptide::data::requireWritableRoot(tmp.path(), 0, 0, true));
}

TEST_CASE("insufficient free space is an error", "[data_paths]") {
    const TempDir tmp;
    const auto huge = std::numeric_limits<std::uintmax_t>::max() / 2;
    CHECK_THROWS_AS(riptide::data::requireWritableRoot(tmp.path(), huge, 0, true), DataPathError);
    CHECK_THROWS_AS(riptide::data::requireWritableRoot(tmp.path(), 1, huge * 2, true),
                    DataPathError);
}

TEST_CASE("run size estimate", "[data_paths]") {
    using riptide::data::estimateRunBytes;
    CHECK(estimateRunBytes(1'000'000, 200) == 200'000'000U);
    CHECK(estimateRunBytes(0, 200) == 0U);
    CHECK_THROWS_AS(estimateRunBytes(std::numeric_limits<std::uint64_t>::max(), 2), DataPathError);
}

TEST_CASE("ids are validated before building paths", "[data_paths]") {
    using riptide::data::isValidId;
    CHECK(isValidId("run_2026-11-16.a"));
    CHECK_FALSE(isValidId(""));
    CHECK_FALSE(isValidId(".."));
    CHECK_FALSE(isValidId(".hidden"));
    CHECK_FALSE(isValidId("-rf"));
    CHECK_FALSE(isValidId("a/b"));
    CHECK_FALSE(isValidId("a b"));
    CHECK_FALSE(isValidId(std::string(riptide::data::kMaxIdLength + 1, 'a')));

    const fs::path root{"/data"};
    CHECK(riptide::data::simRunDir(root, "r1") == fs::path{"/data/sim/r1"});
    CHECK(riptide::data::analysisDir(root, "r1", "a1") == fs::path{"/data/analysis/r1/a1"});
    CHECK_THROWS_AS(riptide::data::simRunDir(root, "../etc"), DataPathError);
    CHECK_THROWS_AS(riptide::data::analysisDir(root, "r1", "x/y"), DataPathError);
    CHECK(riptide::data::validationDir(root, "r1", "c1") == fs::path{"/data/validation/r1/c1"});
    CHECK_THROWS_AS(riptide::data::validationDir(root, "..", "c1"), DataPathError);
    CHECK_THROWS_AS(riptide::data::validationDir(root, "r1", ".c"), DataPathError);
    CHECK(riptide::data::statsDir(root, "cal1") == fs::path{"/data/stats/cal1"});
    CHECK_THROWS_AS(riptide::data::statsDir(root, "../x"), DataPathError);
    CHECK_THROWS_AS(riptide::data::statsDir(root, ""), DataPathError);
}
