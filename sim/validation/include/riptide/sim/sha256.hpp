#pragma once

// SHA-256 (FIPS 180-4) of byte streams and files. The physics checks record
// the digest of every evaluated file they read, so a result can always be
// traced to the exact reference data; computing it in-process avoids
// depending on an external tool at run time.

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <istream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace riptide::sim::checks {

class Sha256Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Incremental hasher: update() any number of times, then digest() once.
class Sha256 {
public:
    static constexpr std::size_t kBlockBytes = 64;
    static constexpr std::size_t kDigestBytes = 32;

    Sha256();
    void update(const unsigned char* data, std::size_t size);
    void update(std::string_view bytes);
    // Pads the message and returns the digest; the hasher must not be reused.
    [[nodiscard]] std::array<std::uint8_t, kDigestBytes> digest();

private:
    void compress(const unsigned char* block);

    std::array<std::uint32_t, 8> state_{};
    std::array<unsigned char, kBlockBytes> buffer_{};
    std::size_t buffered_ = 0;
    std::uint64_t total_bytes_ = 0;
    bool finished_ = false;
};

// Lower-case hexadecimal form of a digest, as printed by sha256sum.
[[nodiscard]] std::string toHex(const std::array<std::uint8_t, Sha256::kDigestBytes>& digest);

// Digest of everything left in the stream, read in blocks.
[[nodiscard]] std::string sha256Hex(std::istream& input);

// Digest of a file; throws Sha256Error if it cannot be opened or read.
[[nodiscard]] std::string sha256File(const std::filesystem::path& path);

} // namespace riptide::sim::checks
