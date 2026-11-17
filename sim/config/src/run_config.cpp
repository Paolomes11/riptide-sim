#include "riptide/sim/run_config.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <string>
#include <string_view>

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

SourceConfig parseSource(const json& obj) {
    constexpr std::string_view where = "source";
    if (!obj.is_object()) {
        fail(where, "expected a JSON object");
    }
    const auto& mode = member(obj, where, "mode");
    if (!mode.is_string()) {
        fail("source.mode", "expected a string");
    }
    SourceConfig source;
    if (mode.get<std::string>() == "mono") {
        requireOnlyKeys(obj, where, {"mode", "En_MeV", "direction", "position_mm"});
        source.mode = SourceMode::Mono;
        source.En_min = finiteNumber(member(obj, where, "En_MeV"), "source.En_MeV") * units::MeV;
        source.En_max = source.En_min;
        checkEnergy(source.En_min, "source.En_MeV");
    } else if (mode.get<std::string>() == "flat") {
        requireOnlyKeys(obj, where,
                        {"mode", "En_min_MeV", "En_max_MeV", "direction", "position_mm"});
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

    const auto direction = vec3(member(obj, where, "direction"), "source.direction", 1.0);
    const double length = norm(direction);
    if (!(std::abs(length - 1.0) <= kDirectionNormTolerance)) {
        fail("source.direction", "expected a unit vector");
    }
    source.direction = {direction.x / length, direction.y / length, direction.z / length};
    source.position = vec3(member(obj, where, "position_mm"), "source.position_mm", units::mm);
    return source;
}

} // namespace

RunConfig parseRunConfig(const json& doc) {
    constexpr std::string_view where = "config";
    requireOnlyKeys(doc, where, {"config_version", "geometry", "source", "physics_list"});
    RunConfig config;

    const auto& version = member(doc, where, "config_version");
    if (!version.is_number_integer() || version.get<int>() != kConfigVersion) {
        fail("config_version", "expected " + std::to_string(kConfigVersion));
    }
    config.config_version = kConfigVersion;
    config.geometry = parseGeometry(member(doc, where, "geometry"));
    config.source = parseSource(member(doc, where, "source"));
    // The gun must start outside the scintillator, otherwise the "first
    // interaction inside the scintillator" of the truth record is ill-defined.
    if (insideBox(config.source.position, config.geometry.scint_size)) {
        fail("source.position_mm", "the source must lie outside the scintillator");
    }

    const auto& physics = member(doc, where, "physics_list");
    if (!physics.is_string() || physics.get<std::string>().empty()) {
        fail("physics_list", "expected a non-empty string");
    }
    config.physics_list = physics.get<std::string>();
    return config;
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
    return parseRunConfig(doc);
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
    source["direction"] = array(s.direction, 1.0);
    source["position_mm"] = array(s.position, units::mm);
    return {{"config_version", config.config_version},
            {"geometry", {{"scint_size_mm", array(config.geometry.scint_size, units::mm)}}},
            {"source", source},
            {"physics_list", config.physics_list}};
}

} // namespace riptide::sim
