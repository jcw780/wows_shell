#ifndef _CONTROL_INDICES_
#define _CONTROL_INDICES_

namespace shell{
namespace impact {
static constexpr unsigned int maxColumns = 13;
static constexpr unsigned int maxColumnsFit = 7;
enum impactDataIndex {
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
static_assert(effectivePenetrationDeckNormalized == (maxColumns - 1),
              "Invalid standard columns");
}  // namespace impact

namespace angle {
static constexpr unsigned int maxColumns = 9;
enum angleDataIndex {
    distance,
    ricochetAngle0Radians,
    ricochetAngle0Degrees,
    ricochetAngle1Radians,
    ricochetAngle1Degrees,
    armorRadians,
    armorDegrees,
    fuseRadians,
    fuseDegrees
};
static_assert(fuseDegrees == (maxColumns - 1), "Invalid angle columns");
}  // namespace angle

namespace post {
static constexpr unsigned int maxColumns = 6;
enum postPenDataIndex { angle, distance, x, y, z, xwf };
static_assert(xwf == (maxColumns - 1), "Invaild postpen columns");
}  // namespace post

enum numerical { forwardEuler, rungeKutta2, rungeKutta4, adamsBashforth5 };
}

#endif