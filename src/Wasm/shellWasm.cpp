#include <emscripten/bind.h>

#include <algorithm>
#include <cstddef>
#include <utility>

#include "../shellCPP.hpp"

class shellCombined {
   private:
    shell::shellCalc calc;
    // shell::shell s;
    std::vector<shell::shell> ships;

   public:
    /*shellCombined(const double caliber, const double v0, const double cD,
                  const double mass, const double krupp,
                  const double normalization, const double fuseTime,
                  const double threshold, const double ricochet0,
                  const double ricochet1) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, "ship");
    }*/
    shellCombined(const unsigned int numShips = 1) { resize(numShips); }
    void resize(const unsigned int numShips) { ships.resize(numShips); }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const double nonAP,
                   const int shipIndex) {
        ships[shipIndex].setValues(caliber, v0, cD, mass, krupp, normalization,
                                   fuseTime, threshold, ricochet0, ricochet1,
                                   nonAP, "");
    }

    void setMax(const double max) { calc.set_max(max); }
    void setMin(const double min) { calc.set_min(min); }
    void setPrecision(const double precision) { calc.set_precision(precision); }
    void setX0(const double x0) { calc.set_x0(x0); }
    void setY0(const double y0) { calc.set_y0(y0); }
    void setDtMin(const double dt) { calc.set_dt_min(dt); }
    void setXf0(const double xf0) { calc.set_xf0(xf0); }
    void setYf0(const double yf0) { calc.set_yf0(yf0); }
    void setDtf(const double dtf) { calc.set_dtf(dtf); }

    // Impact Wrappers
    // Default: Adams Bashforth 5

    template <shell::numerical Numerical>
    void calcImpact() {
        for (auto &s : ships) {
            calc.calculateImpact<Numerical, false>(
                s, false,
                1);  // atomics don't work yet
        }
    }

    // Sizes are the same for both ships
    int impactSize() { return ships[0].impactSize; }
    int impactSizeAligned() { return ships[0].impactSizeAligned; }

    std::vector<double> impactData(const int shipIndex) {
        if (ships[shipIndex].completedImpact) {
            return ships[shipIndex].impactData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getImpactPoint(const int i, const int j, const int shipIndex) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_impact(i, j);
    }

    // Angle Data Wrappers

    void calcAngles(const double thickness, const double inclination) {
        for (auto &s : ships) {
            calc.calculateAngles(thickness, inclination, s, 1);
        }
    }

    std::vector<double> angleData(const int shipIndex) {
        if (ships[shipIndex].completedAngles) {
            return ships[shipIndex].angleData;
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    double getAnglePoint(const int row, const int impact, const int shipIndex) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_angle(row, impact);
    }

    // Post Penetration Wrappers

    void calcPostPen(const double thickness, const double inclination,
                     emscripten::val v, const bool changeDirection,
                     const bool fast) {
        std::vector<double> input =
            std::move(emscripten::vecFromJSArray<double>(v));
        for (auto &s : ships) {
            calc.calculatePostPen(thickness, inclination, s, input,
                                  changeDirection, fast,
                                  1);  // atomics don't work yet
        }
    }

    // Sizes are the same for both ships
    int postPenSize() { return ships[0].postPenSize; }

    std::vector<double> postPenData(const int shipIndex) {
        if (ships[shipIndex].completedPostPen) {
            return ships[shipIndex].postPenData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getPostPenPoint(const int i, const int j, const int k,
                           const int shipIndex) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_postPen(i, j, k);
    }

    // Print Functions

    void printImpact(const int shipIndex) {
        if (ships[shipIndex].completedImpact) {
            ships[shipIndex].printImpactData();
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    void printAngles(const int shipIndex) {
        if (ships[shipIndex].completedAngles) {
            ships[shipIndex].printAngleData();
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    void printPostPen(const int shipIndex) {
        if (ships[shipIndex].completedPostPen) {
            ships[shipIndex].printPostPenData();
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};

// Compile option
// emcc --bind -o shellWasm.wasm.js shellWasm.cpp --std=c++17 -O3 -s
// ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 emcc --bind -o shellWasm.js
// shellWasm.cpp --std=c++17 -O3 -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 -s
// USE_PTHREADS=1 -s WASM_MEM_MAX=100Mb

// Vectorized
// emcc --bind -o shellWasmV.wasm.js shellWasm.cpp --std=c++17 -O3 -s
// ASSERTIONS=1 -s SIMD=1 -msimd128 -s MODULARIZE=1 -s
// 'EXPORT_NAME="ShellWasmV"' -s ENVIRONMENT="web" -s EXPORT_ES6=1 -s
// MALLOC=emmalloc -s ALLOW_MEMORY_GROWTH=1 -s FILESYSTEM=0 -s STRICT=1

// Unvectorized
// emcc --bind -o shellWasm.wasm.js shellWasm.cpp --std=c++17 -O3 -s
// ASSERTIONS=1 -s MODULARIZE=1 -s 'EXPORT_NAME="ShellWasm"' -s
// ENVIRONMENT="web" -s MALLOC=emmalloc -s ALLOW_MEMORY_GROWTH=1 -s FILESYSTEM=0
// -s SINGLE_FILE=1

// Testline: s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
EMSCRIPTEN_BINDINGS(shellWasm) {
    emscripten::class_<shellCombined>("shell")
        .constructor<int>()
        .function("resize", &shellCombined::resize)
        .function("setValues", &shellCombined::setValues)
        .function("setMax", &shellCombined::setMax)
        .function("setMin", &shellCombined::setMin)
        .function("setPrecision", &shellCombined::setPrecision)
        .function("setX0", &shellCombined::setX0)
        .function("setY0", &shellCombined::setY0)
        .function("setDtMin", &shellCombined::setDtMin)
        .function("setXf0", &shellCombined::setXf0)
        .function("setYf0", &shellCombined::setYf0)
        .function("setDtf", &shellCombined::setDtf)

        .function("calcImpact",
                  &shellCombined::calcImpact<shell::numerical::adamsBashforth5>)
        .function("calcImpactAdamsBashforth5",
                  &shellCombined::calcImpact<shell::numerical::adamsBashforth5>)
        .function("calcImpactForwardEuler",
                  &shellCombined::calcImpact<shell::numerical::forwardEuler>)
        .function("calcImpactRungeKutta2",
                  &shellCombined::calcImpact<shell::numerical::rungeKutta2>)
        .function("calcImpactRungeKutta4",
                  &shellCombined::calcImpact<shell::numerical::rungeKutta4>)

        .function("getImpactPoint", &shellCombined::getImpactPoint)
        .function("impactData", &shellCombined::impactData)
        .function("getImpactSize", &shellCombined::impactSize)
        .function("getImpactSizeAligned", &shellCombined::impactSizeAligned)
        .function("calcAngles", &shellCombined::calcAngles)
        .function("angleData", &shellCombined::angleData)
        .function("getAnglePoint", &shellCombined::getAnglePoint)
        .function("calcPostPen", &shellCombined::calcPostPen)
        .function("getPostPenPoint", &shellCombined::getPostPenPoint)
        .function("postPenData", &shellCombined::postPenData)
        .function("getPostPenSize", &shellCombined::postPenSize)
        .function("printImpact", &shellCombined::printImpact)
        .function("printAngles", &shellCombined::printAngles)
        .function("printPostPen", &shellCombined::printPostPen);

    emscripten::register_vector<double>("vector<double>");

    // Enums
    emscripten::enum_<shell::impact::impactIndices>("impactIndices")
        .value("distance", shell::impact::impactIndices::distance)
        .value("launchA", shell::impact::impactIndices::launchAngle)
        .value("impactAHR",
               shell::impact::impactIndices::impactAngleHorizontalRadians)
        .value("impactAHD",
               shell::impact::impactIndices::impactAngleHorizontalDegrees)
        .value("impactV", shell::impact::impactIndices::impactVelocity)
        .value("rawPen", shell::impact::impactIndices::rawPenetration)
        .value("ePenH",
               shell::impact::impactIndices::effectivePenetrationHorizontal)
        .value("ePenHN", shell::impact::impactIndices::
                             effectivePenetrationHorizontalNormalized)
        .value("impactADD",
               shell::impact::impactIndices::impactAngleDeckDegrees)
        .value("ePenD", shell::impact::impactIndices::effectivePenetrationDeck)
        .value("ePenDN",
               shell::impact::impactIndices::effectivePenetrationDeckNormalized)
        .value("tToTarget", shell::impact::impactIndices::timeToTarget)
        .value("tToTargetA",
               shell::impact::impactIndices::timeToTargetAdjusted);

    emscripten::enum_<shell::angle::angleIndices>("angleIndices")
        .value("distance", shell::angle::angleIndices::distance)
        .value("ra0", shell::angle::angleIndices::ricochetAngle0Radians)
        .value("ra0D", shell::angle::angleIndices::ricochetAngle0Degrees)
        .value("ra1", shell::angle::angleIndices::ricochetAngle1Radians)
        .value("ra1D", shell::angle::angleIndices::ricochetAngle1Degrees)
        .value("armor", shell::angle::angleIndices::armorRadians)
        .value("armorD", shell::angle::angleIndices::armorDegrees)
        .value("fuse", shell::angle::angleIndices::fuseRadians)
        .value("fuseD", shell::angle::angleIndices::fuseDegrees);

    emscripten::enum_<shell::post::postPenIndices>("postPenIndices")
        .value("angle", shell::post::postPenIndices::angle)
        .value("distance", shell::post::postPenIndices::distance)
        .value("x", shell::post::postPenIndices::x)
        .value("y", shell::post::postPenIndices::y)
        .value("z", shell::post::postPenIndices::z)
        .value("xwf", shell::post::postPenIndices::xwf);
};