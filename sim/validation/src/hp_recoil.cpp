#include "riptide/sim/hp_recoil.hpp"

#include <cmath>
#include <string>

namespace riptide::sim::checks {

namespace {

constexpr int kMinRepFlag = 1;
constexpr int kMaxRepFlag = 3;

} // namespace

HpElasticHeader parseHpElasticHeader(std::istream& in) {
    HpElasticHeader h;
    if (!(in >> h.rep_flag >> h.awr >> h.frame_flag)) {
        throw HpHeaderError("HP elastic header: expected 'repFlag targetMass frameFlag'");
    }
    if (h.rep_flag < kMinRepFlag || h.rep_flag > kMaxRepFlag) {
        throw HpHeaderError("HP elastic header: unsupported repFlag " + std::to_string(h.rep_flag));
    }
    if (!std::isfinite(h.awr) || !(h.awr > 0.0)) {
        throw HpHeaderError("HP elastic header: target mass must be positive");
    }
    if (h.frame_flag != kHpFrameCm) {
        throw HpHeaderError("HP elastic header: frameFlag " + std::to_string(h.frame_flag) +
                            ", only the CM frame (2) is modelled");
    }
    return h;
}

} // namespace riptide::sim::checks
