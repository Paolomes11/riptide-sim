#include "riptide/sim/run_config.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <string>
#include <string_view>

#include "riptide/sim/source_sampling.hpp"

namespace riptide::sim {

namespace {

using nlohmann::json;

[[noreturn]] void fail(std::string_view where, std::string_view what) {
    throw ConfigError(std::string{where} + ": " + std::string{what});
}

// Unknown keys are errors: a misspelt optional key would otherwise be ignored.
void requireOnlyKeys(const json& obj, std::string_view where,
                     std::initializer_list<std::string_view> allowed) {
    if (!obj.is_object()) {
        fail(where, "expected a JSON object");
    }
    for (const auto& item : obj.items()) {
        bool known = false;
        for (const auto key : allowed) {
            known = known || item.key() == key;
        }
        if (!known) {
            fail(where, "unknown key '" + item.key() + "'");
        }
    }
}

// `key` is a C string, not std::string: a temporary string argument makes GCC
// flag the returned reference with a -Wdangling-reference false positive.
const json& member(const json& obj, std::string_view where, const char* key) {
    const auto it = obj.find(key);
    if (it == obj.end()) {
        fail(where, std::string{"missing key '"} + key + "'");
    }
    return *it;
}

double finiteNumber(const json& value, std::string_view where) {
    if (!value.is_number()) {
        fail(where, "expected a number");
    }
    const auto x = value.get<double>();
    if (!std::isfinite(x)) {
        fail(where, "expected a finite number");
    }
    return x;
}

// Every length and energy goes through the unit constant, so a change of the
// project unit system would not silently reinterpret existing files.
Vec3 vec3(const json& value, std::string_view where, double unit) {
    if (!value.is_array() || value.size() != 3) {
        fail(where, "expected an array of 3 numbers");
    }
    return {finiteNumber(value[0], where) * unit, finiteNumber(value[1], where) * unit,
            finiteNumber(value[2], where) * unit};
}

double norm(const Vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

GeometryConfig parseGeometry(const json& obj) {
    constexpr std::string_view where = "geometry";
    requireOnlyKeys(obj, where, {"scint_size_mm"});
    GeometryConfig geometry;
    geometry.scint_size =
        vec3(member(obj, where, "scint_size_mm"), "geometry.scint_size_mm", units::mm);
    const auto& s = geometry.scint_size;
    if (!(s.x > 0.0 && s.y > 0.0 && s.z > 0.0)) {
        fail("geometry.scint_size_mm", "edge lengths must be positive");
    }
    return geometry;
}

void checkEnergy(double En, std::string_view where) {
    if (!(En > 0.0)) {
        fail(where, "energy must be positive");
    }
    if (En > kNeutronHpMaxEnergy) {
        fail(where, "energy above the 20 MeV limit of the evaluated neutron data");
    }
}

bool insideBox(const Vec3& p, const Vec3& size) {
    return std::abs(p.x) <= 0.5 * size.x && std::abs(p.y) <= 0.5 * size.y &&
           std::abs(p.z) <= 0.5 * size.z;
}

Vec3 unitVector(const json& value, std::string_view where) {
    const auto direction = vec3(value, where, 1.0);
    const double length = norm(direction);
    if (!(std::abs(length - 1.0) <= kDirectionNormTolerance)) {
        fail(where, "expected a unit vector");
    }
    return {direction.x / length, direction.y / length, direction.z / length};
}

// Energy spectrum of the source, common to both beam modes.
void parseSpectrum(const json& obj, SourceConfig& source) {
    constexpr std::string_view where = "source";
    const auto& mode = member(obj, where, "mode");
    if (!mode.is_string()) {
        fail("source.mode", "expected a string");
    }
    if (mode.get<std::string>() == "mono") {
        source.mode = SourceMode::Mono;
        source.En_min = finiteNumber(member(obj, where, "En_MeV"), "source.En_MeV") * units::MeV;
        source.En_max = source.En_min;
        checkEnergy(source.En_min, "source.En_MeV");
    } else if (mode.get<std::string>() == "flat") {
        source.mode = SourceMode::Flat;
        source.En_min =
            finiteNumber(member(obj, where, "En_min_MeV"), "source.En_min_MeV") * units::MeV;
        source.En_max =
            finiteNumber(member(obj, where, "En_max_MeV"), "source.En_max_MeV") * units::MeV;
        checkEnergy(source.En_min, "source.En_min_MeV");
        checkEnergy(source.En_max, "source.En_max_MeV");
        if (!(source.En_min < source.En_max)) {
            fail(where, "En_min_MeV must be smaller than En_max_MeV");
        }
    } else {
        fail("source.mode", "expected \"mono\" or \"flat\"");
    }
}

// Keys of the spectrum for the mode in the file, so that requireOnlyKeys can
// reject keys of the other mode.
bool isFlat(const json& obj) {
    const auto it = obj.find("mode");
    return it != obj.end() && it->is_string() && it->get<std::string>() == "flat";
}

// `keys` is the source object without the beam key, which only version 3 has.
void parsePencil(const json& keys, SourceConfig& source) {
    constexpr std::string_view where = "source";
    if (isFlat(keys)) {
        requireOnlyKeys(keys, where,
                        {"mode", "En_min_MeV", "En_max_MeV", "direction", "position_mm"});
    } else {
        requireOnlyKeys(keys, where, {"mode", "En_MeV", "direction", "position_mm"});
    }
    source.beam = BeamMode::Pencil;
    source.direction = unitVector(member(keys, where, "direction"), "source.direction");
    source.position = vec3(member(keys, where, "position_mm"), "source.position_mm", units::mm);
}

void parseFarField(const json& obj, SourceConfig& source) {
    constexpr std::string_view where = "source";
    if (isFlat(obj)) {
        requireOnlyKeys(obj, where,
                        {"beam", "mode", "En_min_MeV", "En_max_MeV", "direction", "margin_mm"});
    } else {
        requireOnlyKeys(obj, where, {"beam", "mode", "En_MeV", "direction", "margin_mm"});
    }
    source.beam = BeamMode::FarField;
    const auto& direction = member(obj, where, "direction");
    if (direction.is_string()) {
        if (direction.get<std::string>() != "isotropic") {
            fail("source.direction", "expected a unit vector or \"isotropic\"");
        }
        source.isotropic = true;
    } else {
        source.direction = unitVector(direction, "source.direction");
    }
    source.margin = finiteNumber(member(obj, where, "margin_mm"), "source.margin_mm") * units::mm;
    // A positive margin keeps every start point strictly outside the box.
    if (!(source.margin > 0.0)) {
        fail("source.margin_mm", "expected a positive length");
    }
}

SourceConfig parseSource(const json& obj, int config_version) {
    constexpr std::string_view where = "source";
    if (!obj.is_object()) {
        fail(where, "expected a JSON object");
    }
    SourceConfig source;
    parseSpectrum(obj, source);
    if (config_version < kFirstBeamConfigVersion) {
        parsePencil(obj, source);
        return source;
    }
    const auto& beam = member(obj, where, "beam");
    if (!beam.is_string()) {
        fail("source.beam", "expected a string");
    }
    if (beam.get<std::string>() == "pencil") {
        json keys = obj;
        keys.erase("beam");
        parsePencil(keys, source);
    } else if (beam.get<std::string>() == "far_field") {
        parseFarField(obj, source);
    } else {
        fail("source.beam", "expected \"pencil\" or \"far_field\"");
    }
    return source;
}

ScintillatorConfig parseScintillator(const json& obj) {
    constexpr std::string_view where = "scintillator";
    requireOnlyKeys(obj, where, {"birks_kB_mm_per_MeV", "birks_source"});
    ScintillatorConfig scintillator;
    scintillator.birks_kB = finiteNumber(member(obj, where, "birks_kB_mm_per_MeV"),
                                         "scintillator.birks_kB_mm_per_MeV") *
                            kBirksUnit;
    // Zero is allowed on purpose: it switches quenching off, which is the
    // deterministic oracle visible energy == deposited energy.
    if (scintillator.birks_kB < 0.0) {
        fail("scintillator.birks_kB_mm_per_MeV", "the Birks constant must not be negative");
    }
    const auto& source = member(obj, where, "birks_source");
    if (!source.is_string() || source.get<std::string>().empty()) {
        fail("scintillator.birks_source", "expected a non-empty string");
    }
    scintillator.birks_source = source.get<std::string>();
    return scintillator;
}

StepFunction parseStepFunction(const json& obj) {
    constexpr std::string_view where = "proton_step_function";
    requireOnlyKeys(obj, where, {"dRoverRange", "finalRange_mm"});
    StepFunction step;
    step.dRoverRange =
        finiteNumber(member(obj, where, "dRoverRange"), "proton_step_function.dRoverRange");
    step.finalRange =
        finiteNumber(member(obj, where, "finalRange_mm"), "proton_step_function.finalRange_mm") *
        units::mm;
    if (!(step.dRoverRange > 0.0 && step.dRoverRange <= kMaxDRoverRange)) {
        fail("proton_step_function.dRoverRange", "expected a value in (0, 1]");
    }
    if (!(step.finalRange > 0.0)) {
        fail("proton_step_function.finalRange_mm", "expected a positive length");
    }
    return step;
}

int parseVersion(const json& doc) {
    const auto& version = member(doc, "config", "config_version");
    if (!version.is_number_integer() || version.get<int>() < kOldestConfigVersion ||
        version.get<int>() > kConfigVersion) {
        fail("config_version", "expected a version from " + std::to_string(kOldestConfigVersion) +
                                   " to " + std::to_string(kConfigVersion));
    }
    return version.get<int>();
}

} // namespace

RunConfig parseRunConfig(const json& doc) {
    constexpr std::string_view where = "config";
    if (!doc.is_object()) {
        fail(where, "expected a JSON object");
    }
    RunConfig config;
    config.config_version = parseVersion(doc);
    if (config.config_version >= kFirstBirksConfigVersion) {
        requireOnlyKeys(doc, where,
                        {"config_version", "geometry", "source", "physics_list", "scintillator",
                         "proton_step_function"});
        config.scintillator = parseScintillator(member(doc, where, "scintillator"));
        config.proton_step_function = parseStepFunction(member(doc, where, "proton_step_function"));
    } else {
        requireOnlyKeys(doc, where, {"config_version", "geometry", "source", "physics_list"});
    }
    config.geometry = parseGeometry(member(doc, where, "geometry"));
    config.source = parseSource(member(doc, where, "source"), config.config_version);
    // The gun must start outside the scintillator, otherwise the "first
    // interaction inside the scintillator" of the truth record is ill-defined.
    // The far-field disk lies outside by construction (source_sampling.hpp).
    if (config.source.beam == BeamMode::Pencil &&
        insideBox(config.source.position, config.geometry.scint_size)) {
        fail("source.position_mm", "the source must lie outside the scintillator");
    }

    const auto& physics = member(doc, where, "physics_list");
    if (!physics.is_string() || physics.get<std::string>().empty()) {
        fail("physics_list", "expected a non-empty string");
    }
    config.physics_list = physics.get<std::string>();
    return config;
}

Vec3 scintHalfSize(const GeometryConfig& geometry) {
    const auto& s = geometry.scint_size;
    return {0.5 * s.x, 0.5 * s.y, 0.5 * s.z};
}

double worldHalfSize(const RunConfig& config) {
    const Vec3 half = scintHalfSize(config.geometry);
    const double extent = std::max({half.x, half.y, half.z, sourceExtent(config.source, half)});
    return kWorldScale * extent + kWorldMargin;
}

RunConfig loadRunConfig(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw ConfigError("cannot open config file: " + path.string());
    }
    json doc;
    try {
        in >> doc;
    } catch (const json::parse_error& e) {
        throw ConfigError("invalid JSON in " + path.string() + ": " + e.what());
    }
    auto config = parseRunConfig(doc);
    if (config.config_version != kConfigVersion) {
        throw ConfigError(path.string() + ": config_version " +
                          std::to_string(config.config_version) +
                          " is only readable, a new run needs " + std::to_string(kConfigVersion));
    }
    return config;
}

json toJson(const RunConfig& config) {
    const auto array = [](const Vec3& v, double unit) {
        return json::array({v.x / unit, v.y / unit, v.z / unit});
    };
    json source;
    const auto& s = config.source;
    if (s.mode == SourceMode::Mono) {
        source = {{"mode", "mono"}, {"En_MeV", s.En_min / units::MeV}};
    } else {
        source = {{"mode", "flat"},
                  {"En_min_MeV", s.En_min / units::MeV},
                  {"En_max_MeV", s.En_max / units::MeV}};
    }
    if (config.config_version >= kFirstBeamConfigVersion) {
        source["beam"] = s.beam == BeamMode::Pencil ? "pencil" : "far_field";
    }
    if (s.beam == BeamMode::Pencil) {
        source["direction"] = array(s.direction, 1.0);
        source["position_mm"] = array(s.position, units::mm);
    } else {
        source["direction"] = s.isotropic ? json("isotropic") : array(s.direction, 1.0);
        source["margin_mm"] = s.margin / units::mm;
    }
    json doc = {{"config_version", config.config_version},
                {"geometry", {{"scint_size_mm", array(config.geometry.scint_size, units::mm)}}},
                {"source", source},
                {"physics_list", config.physics_list}};
    if (config.scintillator) {
        doc["scintillator"] = {{"birks_kB_mm_per_MeV", config.scintillator->birks_kB / kBirksUnit},
                               {"birks_source", config.scintillator->birks_source}};
    }
    if (config.proton_step_function) {
        doc["proton_step_function"] = {
            {"dRoverRange", config.proton_step_function->dRoverRange},
            {"finalRange_mm", config.proton_step_function->finalRange / units::mm}};
    }
    return doc;
}

} // namespace riptide::sim
