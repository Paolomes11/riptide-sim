// SHA-256 against the FIPS 180-2 example messages and padding edge cases.
// Expected digests were produced with GNU coreutils sha256sum.

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "riptide/sim/sha256.hpp"

namespace checks = riptide::sim::checks;

namespace {

std::string hashOf(const std::string& message) {
    checks::Sha256 hasher;
    hasher.update(message);
    return checks::toHex(hasher.digest());
}

} // namespace

TEST_CASE("sha256 of the FIPS example messages", "[sha256]") {
    CHECK(hashOf("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    CHECK(hashOf("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    CHECK(hashOf("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
          "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    constexpr std::size_t kMillion = 1000000;
    CHECK(hashOf(std::string(kMillion, 'a')) ==
          "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

TEST_CASE("sha256 padding at the block boundaries", "[sha256]") {
    // 55 bytes: the length still fits in the same block; 64: one full block.
    CHECK(hashOf(std::string(55, 'x')) ==
          "d5e285683cd4efc02d021a5c62014694958901005d6f71e89e0989fac77e4072");
    CHECK(hashOf(std::string(64, 'x')) ==
          "7ce100971f64e7001e8fe5a51973ecdfe1ced42befe7ee8d5fd6219506b5393c");
}

TEST_CASE("sha256 does not depend on how the input is split", "[sha256]") {
    const std::string message = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    checks::Sha256 hasher;
    for (const char c : message) {
        hasher.update(std::string(1, c));
    }
    CHECK(checks::toHex(hasher.digest()) == hashOf(message));

    std::istringstream stream(message);
    CHECK(checks::sha256Hex(stream) == hashOf(message));
}

TEST_CASE("sha256 refuses reuse and missing files", "[sha256]") {
    checks::Sha256 hasher;
    (void)hasher.digest();
    CHECK_THROWS_AS(hasher.update("a"), checks::Sha256Error);
    CHECK_THROWS_AS(hasher.digest(), checks::Sha256Error);
    CHECK_THROWS_AS(checks::sha256File("/nonexistent/riptide/sha256/file"), checks::Sha256Error);
}

TEST_CASE("sha256 of a file", "[sha256]") {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "riptide_sha256_test.bin";
    {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out << "abc";
    }
    CHECK(checks::sha256File(path) ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    std::filesystem::remove(path);
}
