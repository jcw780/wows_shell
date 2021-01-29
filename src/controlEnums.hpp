#pragma once
#include <cstdint>
#include <cstdlib>
#include <type_traits>

namespace wows_shell {

// https://stackoverflow.com/questions/8357240/how-to-automatically-convert-strongly-typed-enum-into-int
template <typename E>
constexpr typename std::underlying_type<E>::type toUnderlying(E e) noexcept {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

namespace impact {
static constexpr std::size_t maxColumns = 13;
static constexpr std::size_t maxColumnsFit = 7;
enum class impactIndices {
    distance,
    launchAngle,
    impactAngleHorizontalRadians,  // Negative for Falling
    impactAngleHorizontalDegrees,  // Positive for Falling
    impactVelocity,
    timeToTarget,
    timeToTargetAdjusted,  // Adjusted for in game shell time
    rawPenetration,
    effectivePenetrationHorizontal,
    effectivePenetrationHorizontalNormalized,
    impactAngleDeckDegrees,
    effectivePenetrationDeck,
    effectivePenetrationDeckNormalized,
};
using indexT = typename std::underlying_type<impactIndices>::type;
static_assert(toUnderlying(impactIndices::effectivePenetrationDeckNormalized) ==
                  (maxColumns - 1),
              "Invalid standard columns");
}  // namespace impact

namespace angle {
static constexpr std::size_t maxColumns = 8;
enum class angleIndices {
    ricochetAngle0Radians,
    ricochetAngle0Degrees,
    ricochetAngle1Radians,
    ricochetAngle1Degrees,
    armorRadians,
    armorDegrees,
    fuseRadians,
    fuseDegrees
};
using indexT = typename std::underlying_type<angleIndices>::type;
static_assert(toUnderlying(angleIndices::fuseDegrees) == (maxColumns - 1),
              "Invalid angle columns");
}  // namespace angle

namespace dispersion {
static constexpr std::size_t maxColumns = 9;
enum class dispersionIndices {
    maxHorizontal,
    standardHorizontal,
    halfHorizontal,
    maxVertical,
    standardVertical,
    halfVertical,
    maxArea,
    standardArea,
    halfArea,
};
static_assert(toUnderlying(dispersionIndices::halfArea) == (maxColumns - 1),
              "Invalid dispersion columns");
enum class verticalTypes { horizontal, normal, vertical };
static_assert(toUnderlying(verticalTypes::vertical) == 2,
              "Invalid vertical types");
};  // namespace dispersion

namespace post {
static constexpr std::size_t maxColumns = 5;
enum class postPenIndices { angle, x, y, z, xwf };
using indexT = typename std::underlying_type<postPenIndices>::type;
static_assert(toUnderlying(postPenIndices::xwf) == (maxColumns - 1),
              "Invaild postpen columns");
}  // namespace post

namespace calculateType {
enum class calcIndices { impact, angle, dispersion, post };
static_assert(toUnderlying(calcIndices::post) == 3, "Invalid data indices");
}  // namespace calculateType

enum class numerical {
    forwardEuler,
    rungeKutta2,
    rungeKutta4,
    adamsBashforth5
};

template <numerical Numerical>
static constexpr bool isMultistep() {
    if constexpr (Numerical == numerical::adamsBashforth5) {
        return true;
    } else {
        return false;
    }
}
}  // namespace wows_shell
