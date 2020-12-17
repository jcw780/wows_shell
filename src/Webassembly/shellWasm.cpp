#include <emscripten/bind.h>

#include <algorithm>
#include <cstddef>
#include <utility>

#include "../shellCPP.hpp"

namespace wows_shell {
namespace pointArray {
// These functions are for producing arrays suitable for chart.js scatter plots
emscripten::val getImpactPointArray(shell &s, const std::size_t xIndex,
                                    const std::size_t yIndex) {
    if (s.completedImpact) {
        emscripten::val points = emscripten::val::array();
        for (std::size_t i = 0; i < s.impactSize; i++) {
            emscripten::val point = emscripten::val::object();
            point.set("x", s.get_impact(i, xIndex));
            point.set("y", s.get_impact(i, yIndex));
            points.call<void>("push", point);
        }
        return points;
    } else {
        throw std::runtime_error("Impact data not generated");
    }
}

emscripten::val getAnglePointArray(shell &s, const std::size_t xIndex,
                                   const std::size_t yIndex) {
    if (s.completedAngles) {
        emscripten::val points = emscripten::val::array();
        for (std::size_t i = 0; i < s.impactSize; i++) {
            emscripten::val point = emscripten::val::object();
            point.set("x", s.get_angle(i, xIndex));
            point.set("y", s.get_angle(i, yIndex));
            points.call<void>("push", point);
        }
        return points;
    } else {
        throw std::runtime_error("Impact data not generated");
    }
}

emscripten::val getPostPenPointArray(shell &s, const std::size_t angle,
                                     const std::size_t xIndex,
                                     const std::size_t yIndex) {
    if (s.completedPostPen) {
        emscripten::val points = emscripten::val::array();
        for (std::size_t i = 0; i < s.impactSize; i++) {
            emscripten::val point = emscripten::val::object();
            point.set("x", s.get_postPen(i, xIndex, angle));
            point.set("y", s.get_postPen(i, yIndex, angle));
            points.call<void>("push", point);
        }
        return points;
    } else {
        throw std::runtime_error("Impact data not generated");
    }
}

emscripten::val getPostPenPointArrayFuseStatus(shell &s,
                                               const bool addCondition,
                                               const std::size_t angle,
                                               const std::size_t xIndex,
                                               const std::size_t yIndex) {
    if (s.completedImpact) {
        emscripten::val points = emscripten::val::array();
        if (addCondition) {
            for (std::size_t i = 0; i < s.impactSize; i++) {
                if (s.get_postPen(i, post::postPenIndices::xwf, angle) >= 0) {
                    emscripten::val point = emscripten::val::object();
                    point.set("x", s.get_postPen(i, xIndex, angle));
                    point.set("y", s.get_postPen(i, yIndex, angle));
                    points.call<void>("push", point);
                }
            }
        } else {
            for (std::size_t i = 0; i < s.impactSize; i++) {
                if (s.get_postPen(i, post::postPenIndices::xwf, angle) < 0) {
                    emscripten::val point = emscripten::val::object();
                    point.set("x", s.get_postPen(i, xIndex, angle));
                    point.set("y", s.get_postPen(i, yIndex, angle));
                    points.call<void>("push", point);
                }
            }
        }
        return points;
    } else {
        throw std::runtime_error("Post pen data not generated");
    }
}

emscripten::val getImpactSizedPointArray(shell &s, emscripten::val params1,
                                         emscripten::val params2) {
    emscripten::val points = emscripten::val::array();
    double *startX, *startY;
    auto setup = [&](emscripten::val &param, double *&tgt) {
        std::size_t index = param[0].as<std::size_t>();
        switch (index) {
            case toUnderlying(calculateType::calcIndices::impact):
                if (s.completedImpact) {
                    std::size_t index = param[1].as<std::size_t>();
                    if (0 <= index && index < impact::maxColumns) {
                        tgt = s.get_impactPtr(0, index);
                    } else {
                        throw std::runtime_error("Invalid impact index");
                    }
                } else {
                    throw std::runtime_error("Impact data not generated");
                }
                break;
            case toUnderlying(calculateType::calcIndices::angle):
                if (s.completedAngles) {
                    std::size_t index = param[1].as<std::size_t>();
                    if (0 <= index && index < angle::maxColumns) {
                        tgt = s.get_anglePtr(0, index);
                    } else {
                        throw std::runtime_error("Invalid angle index");
                    }
                } else {
                    throw std::runtime_error("Angle data not generated");
                }
                break;
            case toUnderlying(calculateType::calcIndices::dispersion):
                if (s.completedDispersion) {
                    std::size_t index = param[1].as<std::size_t>();
                    if (0 <= index && index < dispersion::maxColumns) {
                        tgt = s.get_dispersionPtr(0, index);
                    } else {
                        throw std::runtime_error("Invalid dispersion index");
                    }
                } else {
                    throw std::runtime_error("Dispersion data not generated");
                }
                break;
            case toUnderlying(calculateType::calcIndices::post):
                if (s.completedPostPen) {
                    std::size_t index = param[1].as<std::size_t>();
                    if (0 <= index && index < post::maxColumns) {
                        tgt = s.get_postPenPtr(0, index,
                                               param[2].as<std::size_t>());
                    } else {
                        throw std::runtime_error("Invalid post index");
                    }
                } else {
                    throw std::runtime_error("Post data not generated");
                }
                break;
            default:
                throw std::runtime_error("Invalid selection");
        }
    };
    setup(params1, startX);
    setup(params2, startY);

    for (std::size_t i = 0; i < s.impactSize; ++i) {
        emscripten::val point = emscripten::val::object();
        point.set("x", *(startX + i));
        point.set("y", *(startY + i));
        points.call<void>("push", point);
    }

    return points;
}

}  // namespace pointArray

#define ENABLE_SPLIT_SHELL
#ifdef ENABLE_SPLIT_SHELL
class shellWasm {
   public:
    shell s;
    shellWasm(const double caliber, const double v0, const double cD,
              const double mass, const double krupp, const double normalization,
              const double fuseTime, const double threshold,
              const double ricochet0, const double ricochet1,
              const double nonAP, const std::string &name) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, nonAP, name);
    }
    shellWasm(const shellParams &sp, const std::string &name) {
        s.setValues(sp, name);
    }
    shellWasm(const shellParams &sp, const dispersionParams &dp,
              const std::string &name) {
        s.setValues(sp, dp, name);
    }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const double nonAP,
                   const std::string &name) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, nonAP, name);
    }
    void setValues(const shellParams &sp, const std::string &name) {
        s.setValues(sp, name);
    }
    void setValues(const shellParams &sp, const dispersionParams &dp,
                   const std::string &name) {
        s.setValues(sp, dp, name);
    }

    std::size_t impactSize() { return s.impactSize; }
    std::size_t impactSizeAligned() { return s.impactSizeAligned; }

    std::vector<double> impactData() {
        if (s.completedImpact) {
            return s.impactData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    // These functions are for producing arrays suitable for chart.js scatter
    // plots
    emscripten::val getImpactPointArray(const std::size_t xIndex,
                                        const std::size_t yIndex) {
        return pointArray::getImpactPointArray(s, xIndex, yIndex);
    }

    double getImpactPoint(const std::size_t i, const std::size_t j) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return s.get_impact(i, j);
    }

    std::vector<double> angleData() {
        if (s.completedAngles) {
            return s.angleData;
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    double getAnglePoint(const std::size_t row, const std::size_t impact) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return s.get_angle(row, impact);
    }

    // These functions are for producing arrays suitable for chart.js scatter
    // plots
    emscripten::val getAnglePointArray(const std::size_t xIndex,
                                       const std::size_t yIndex) {
        return pointArray::getAnglePointArray(s, xIndex, yIndex);
    }

    std::size_t postPenSize() { return s.postPenSize; }

    std::vector<double> postPenData() {
        if (s.completedPostPen) {
            return s.postPenData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getPostPenPoint(const std::size_t i, const std::size_t j,
                           const std::size_t k) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return s.get_postPen(i, j, k);
    }

    // These functions are for producing arrays suitable for chart.js scatter
    // plots
    emscripten::val getPostPenPointArray(const std::size_t angle,
                                         const std::size_t xIndex,
                                         const std::size_t yIndex) {
        return pointArray::getPostPenPointArray(s, angle, xIndex, yIndex);
    }

    // These functions are for producing arrays suitable for chart.js scatter
    // plots
    emscripten::val getPostPenPointArrayFuseStatus(const bool addCondition,
                                                   const std::size_t angle,
                                                   const std::size_t xIndex,
                                                   const std::size_t yIndex) {
        return pointArray::getPostPenPointArrayFuseStatus(
            s, addCondition, angle, xIndex, yIndex);
    }

    void printImpact() {
        if (s.completedImpact) {
            s.printImpactData();
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    void printAngles() {
        if (s.completedAngles) {
            s.printAngleData();
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    void printPostPen() {
        if (s.completedPostPen) {
            s.printPostPenData();
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};

std::string generateShellWasmHash(const shellWasm &s) {
    return generateHash(s.s);
}

emscripten::val getImpactSizedPointArray(shellWasm &s, emscripten::val params1,
                                         emscripten::val params2) {
    return pointArray::getImpactSizedPointArray(s.s, params1, params2);
}

class shellCalcWasm : public shellCalc {
   public:
    shellCalcWasm() = default;

    void setMax(const double max) { set_max(max); }
    void setMin(const double min) { set_min(min); }
    void setPrecision(const double precision) { set_precision(precision); }
    void setX0(const double x0) { set_x0(x0); }
    void setY0(const double y0) { set_y0(y0); }
    void setDtMin(const double dt) { set_dt_min(dt); }
    void setXf0(const double xf0) { set_xf0(xf0); }
    void setYf0(const double yf0) { set_yf0(yf0); }
    void setDtf(const double dtf) { set_dtf(dtf); }

    template <numerical Numerical>
    void calcImpact(shellWasm &sp) {
#ifdef __EMSCRIPTEN_PTHREADS__
        calculateImpact<false, Numerical, false>(sp.s);
#else
        calculateImpact<false, Numerical, false>(sp.s, 1);
#endif
    }

    void calcAngles(shellWasm &sp, const double thickness,
                    const double inclination) {
#ifdef __EMSCRIPTEN_PTHREADS__
        calculateAngles(thickness, inclination, sp.s);
#else
        calculateAngles(thickness, inclination, sp.s, 1);
#endif
    }

    void calcDispersion(shellWasm &sp) {
#ifdef __EMSCRIPTEN_PTHREADS__
        calculateDispersion(sp.s);
#else
        calculateDispersion(sp.s, 1);
#endif
    }

    void calcPostPen(shellWasm &sp, const double thickness,
                     const double inclination, emscripten::val anglesVal,
                     const bool changeDirection, const bool fast) {
        std::vector<double> input =
            emscripten::convertJSArrayToNumberVector<double>(anglesVal);
#ifdef __EMSCRIPTEN_PTHREADS__
        calculatePostPen(thickness, inclination, sp.s, input, changeDirection,
                         fast);
#else
        calculatePostPen(thickness, inclination, sp.s, input, changeDirection,
                         fast, 1);
#endif
    }
};
#endif

#define ENABLE_SHELL_COMBINED
#ifdef ENABLE_SHELL_COMBINED
class shellCombined {
   private:
    shellCalc calc;
    // shell s;
    std::vector<shell> ships;

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

    template <numerical Numerical>
    void calcImpact() {
        for (auto &s : ships) {
#ifdef __EMSCRIPTEN_PTHREADS__
            calc.calculateImpact<Numerical, false>(s, false);
#else
            calc.calculateImpact<Numerical, false>(s, false, 1);
#endif
        }
    }

    // Sizes are the same for both ships
    int impactSize() { return ships[0].impactSize; }
    int impactSizeAligned() { return ships[0].impactSizeAligned; }

    std::vector<double> impactData(const std::size_t shipIndex) {
        if (ships[shipIndex].completedImpact) {
            return ships[shipIndex].impactData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getImpactPoint(const std::size_t shipIndex, const std::size_t i,
                          const std::size_t j) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_impact(i, j);
    }

    // These functions are for producing arrays suitable for chart.js scatter
    // plots
    emscripten::val getImpactPointArray(const std::size_t shipIndex,
                                        const std::size_t xIndex,
                                        const std::size_t yIndex) {
        return pointArray::getImpactPointArray(ships[shipIndex], xIndex,
                                               yIndex);
    }

    // Angle Data Wrappers

    void calcAngles(const double thickness, const double inclination) {
        for (auto &s : ships) {
#ifdef __EMSCRIPTEN_PTHREADS__
            calc.calculateAngles(thickness, inclination, s);
#else
            calc.calculateAngles(thickness, inclination, s, 1);
#endif
        }
    }

    std::vector<double> angleData(const std::size_t shipIndex) {
        if (ships[shipIndex].completedAngles) {
            return ships[shipIndex].angleData;
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    double getAnglePoint(const std::size_t shipIndex, const std::size_t row,
                         const std::size_t impact) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_angle(row, impact);
    }

    emscripten::val getAnglePointArray(const std::size_t shipIndex,
                                       const std::size_t xIndex,
                                       const std::size_t yIndex) {
        return pointArray::getAnglePointArray(ships[shipIndex], xIndex, yIndex);
    }

    // Post Penetration Wrappers

    void calcPostPen(const double thickness, const double inclination,
                     emscripten::val v, const bool changeDirection,
                     const bool fast) {
        std::vector<double> input =
            emscripten::convertJSArrayToNumberVector<double>(v);
        for (auto &s : ships) {
#ifdef __EMSCRIPTEN_PTHREADS__
            calc.calculatePostPen(thickness, inclination, s, input,
                                  changeDirection, fast);
#else
            calc.calculatePostPen(thickness, inclination, s, input,
                                  changeDirection, fast, 1);
#endif
        }
    }

    // Sizes are the same for both ships
    int postPenSize() { return ships[0].postPenSize; }

    std::vector<double> postPenData(const std::size_t shipIndex) {
        if (ships[shipIndex].completedPostPen) {
            return ships[shipIndex].postPenData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getPostPenPoint(const std::size_t shipIndex, const std::size_t i,
                           const std::size_t j, const std::size_t k) {
        // NOT SAFE - PLEASE MAKE SURE YOU ARE NOT OVERFLOWING
        return ships[shipIndex].get_postPen(i, j, k);
    }

    emscripten::val getPostPenPointArray(const std::size_t shipIndex,
                                         const std::size_t angle,
                                         const std::size_t xIndex,
                                         const std::size_t yIndex) {
        return pointArray::getPostPenPointArray(ships[shipIndex], angle, xIndex,
                                                yIndex);
    }

    emscripten::val getPostPenPointArrayFuseStatus(const std::size_t shipIndex,
                                                   const bool addCondition,
                                                   const std::size_t angle,
                                                   const std::size_t xIndex,
                                                   const std::size_t yIndex) {
        return pointArray::getPostPenPointArrayFuseStatus(
            ships[shipIndex], addCondition, angle, xIndex, yIndex);
    }

    // Print Functions

    void printImpact(const std::size_t shipIndex) {
        if (ships[shipIndex].completedImpact) {
            ships[shipIndex].printImpactData();
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    void printAngles(const std::size_t shipIndex) {
        if (ships[shipIndex].completedAngles) {
            ships[shipIndex].printAngleData();
        } else {
            throw std::runtime_error("Angle data not generated");
        }
    }

    void printPostPen(const std::size_t shipIndex) {
        if (ships[shipIndex].completedPostPen) {
            ships[shipIndex].printPostPenData();
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};
#endif

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
    emscripten::value_object<shellParams>("shellParams")
        .field("caliber", &shellParams::caliber)
        .field("v0", &shellParams::v0)
        .field("cD", &shellParams::cD)
        .field("mass", &shellParams::mass)
        .field("krupp", &shellParams::krupp)
        .field("normalization", &shellParams::normalization)
        .field("fuseTime", &shellParams::fuseTime)
        .field("threshold", &shellParams::threshold)
        .field("ricochet0", &shellParams::ricochet0)
        .field("ricochet1", &shellParams::ricochet1)
        .field("nonAP", &shellParams::nonAP);

    emscripten::value_object<dispersionParams>("dispersionParams")
        .field("idealRadius", &dispersionParams::idealRadius)
        .field("minRadius", &dispersionParams::minRadius)
        .field("idealDistance", &dispersionParams::idealDistance)
        .field("taperDistance", &dispersionParams::taperDistance)
        .field("delim", &dispersionParams::delim)
        .field("zeroRadius", &dispersionParams::zeroRadius)
        .field("delimRadius", &dispersionParams::delimRadius)
        .field("maxRadius", &dispersionParams::maxRadius)
        .field("maxDistance", &dispersionParams::maxDistance)
        .field("sigma", &dispersionParams::sigma);

#ifdef ENABLE_SPLIT_SHELL
    emscripten::class_<shellWasm>("shell")
        .constructor<double, double, double, double, double, double, double,
                     double, double, double, double, std::string>()
        .constructor<shellParams, std::string>()
        .constructor<shellParams, dispersionParams, std::string>()
        .function(
            "setValues",
            static_cast<void (shellWasm::*)(
                const double, const double, const double, const double,
                const double, const double, const double, const double,
                const double, const double, const double, const std::string &)>(
                &shellWasm::setValues))
        .function("setValues", static_cast<void (shellWasm::*)(
                                   const shellParams &, const std::string &)>(
                                   &shellWasm::setValues))
        .function("setValues",
                  static_cast<void (shellWasm::*)(
                      const shellParams &, const dispersionParams &,
                      const std::string &)>(&shellWasm::setValues))
        .function("getImpactPoint", &shellWasm::getImpactPoint)
        .function("getImpactPointArray", &shellWasm::getImpactPointArray)
        .function("impactData", &shellWasm::impactData)
        .function("getImpactSize", &shellWasm::impactSize)
        .function("getImpactSizeAligned", &shellWasm::impactSizeAligned)
        .function("angleData", &shellWasm::angleData)
        .function("getAnglePoint", &shellWasm::getAnglePoint)
        .function("getAnglePointArray", &shellWasm::getAnglePointArray)
        .function("postPenData", &shellWasm::postPenData)
        .function("getPostPenPoint", &shellWasm::getPostPenPoint)
        .function("getPostPenPointArray", &shellWasm::getPostPenPointArray)
        .function("getPostPenPointArrayFuseStatus",
                  &shellWasm::getPostPenPointArrayFuseStatus)
        .function("getPostPenSize", &shellWasm::postPenSize)
        .function("printImpact", &shellWasm::printImpact)
        .function("printAngles", &shellWasm::printAngles)
        .function("printPostPen", &shellWasm::printPostPen);

    emscripten::function("generateHash", &generateShellParamHash);
    emscripten::function("generateShellHash", &generateShellWasmHash);
    emscripten::function("getImpactSizedPointArray", &getImpactSizedPointArray);

    emscripten::class_<shellCalcWasm>("shellCalc")
        .constructor()
        .function("setMax", &shellCalcWasm::setMax)
        .function("setMin", &shellCalcWasm::setMin)
        .function("setPrecision", &shellCalcWasm::setPrecision)
        .function("setX0", &shellCalcWasm::setX0)
        .function("setY0", &shellCalcWasm::setY0)
        .function("setDtMin", &shellCalcWasm::setDtMin)
        .function("setXf0", &shellCalcWasm::setXf0)
        .function("setYf0", &shellCalcWasm::setYf0)
        .function("setDtf", &shellCalcWasm::setDtf)
        .function("calcImpact",
                  &shellCalcWasm::calcImpact<numerical::forwardEuler>)
        .function("calcImpactAdamsBashforth5",
                  &shellCalcWasm::calcImpact<numerical::adamsBashforth5>)
        .function("calcImpactForwardEuler",
                  &shellCalcWasm::calcImpact<numerical::forwardEuler>)
        .function("calcImpactRungeKutta2",
                  &shellCalcWasm::calcImpact<numerical::rungeKutta2>)
        .function("calcImpactRungeKutta4",
                  &shellCalcWasm::calcImpact<numerical::rungeKutta4>)
        .function("calcAngles", &shellCalcWasm::calcAngles)
        .function("calcDispersion", &shellCalcWasm::calcDispersion)
        .function("calcPostPen", &shellCalcWasm::calcPostPen);
#endif
#ifdef ENABLE_SHELL_COMBINED
    emscripten::class_<shellCombined>("shellCombined")
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
                  &shellCombined::calcImpact<numerical::adamsBashforth5>)
        .function("calcImpactAdamsBashforth5",
                  &shellCombined::calcImpact<numerical::adamsBashforth5>)
        .function("calcImpactForwardEuler",
                  &shellCombined::calcImpact<numerical::forwardEuler>)
        .function("calcImpactRungeKutta2",
                  &shellCombined::calcImpact<numerical::rungeKutta2>)
        .function("calcImpactRungeKutta4",
                  &shellCombined::calcImpact<numerical::rungeKutta4>)

        .function("getImpactPoint", &shellCombined::getImpactPoint)
        .function("impactData", &shellCombined::impactData)
        .function("getImpactPointArray", &shellCombined::getImpactPointArray)
        .function("getImpactSize", &shellCombined::impactSize)
        .function("getImpactSizeAligned", &shellCombined::impactSizeAligned)

        .function("calcAngles", &shellCombined::calcAngles)
        .function("angleData", &shellCombined::angleData)
        .function("getAnglePoint", &shellCombined::getAnglePoint)
        .function("getAnglePointArray", &shellCombined::getAnglePointArray)

        .function("calcPostPen", &shellCombined::calcPostPen)
        .function("getPostPenPoint", &shellCombined::getPostPenPoint)
        .function("postPenData", &shellCombined::postPenData)
        .function("getPostPenPointArray", &shellCombined::getPostPenPointArray)
        .function("getPostPenPointArrayFuseStatus",
                  &shellCombined::getPostPenPointArrayFuseStatus)
        .function("getPostPenSize", &shellCombined::postPenSize)

        .function("printImpact", &shellCombined::printImpact)
        .function("printAngles", &shellCombined::printAngles)
        .function("printPostPen", &shellCombined::printPostPen);
#endif
    emscripten::register_vector<double>("vector<double>");

    // Enums
    emscripten::enum_<impact::impactIndices>("impactIndices")
        .value("distance", impact::impactIndices::distance)
        .value("launchA", impact::impactIndices::launchAngle)
        .value("impactAHR", impact::impactIndices::impactAngleHorizontalRadians)
        .value("impactAHD", impact::impactIndices::impactAngleHorizontalDegrees)
        .value("impactV", impact::impactIndices::impactVelocity)
        .value("rawPen", impact::impactIndices::rawPenetration)
        .value("ePenH", impact::impactIndices::effectivePenetrationHorizontal)
        .value("ePenHN",
               impact::impactIndices::effectivePenetrationHorizontalNormalized)
        .value("impactADD", impact::impactIndices::impactAngleDeckDegrees)
        .value("ePenD", impact::impactIndices::effectivePenetrationDeck)
        .value("ePenDN",
               impact::impactIndices::effectivePenetrationDeckNormalized)
        .value("tToTarget", impact::impactIndices::timeToTarget)
        .value("tToTargetA", impact::impactIndices::timeToTargetAdjusted);

    emscripten::enum_<angle::angleIndices>("angleIndices")
        .value("distance", angle::angleIndices::distance)
        .value("ra0", angle::angleIndices::ricochetAngle0Radians)
        .value("ra0D", angle::angleIndices::ricochetAngle0Degrees)
        .value("ra1", angle::angleIndices::ricochetAngle1Radians)
        .value("ra1D", angle::angleIndices::ricochetAngle1Degrees)
        .value("armor", angle::angleIndices::armorRadians)
        .value("armorD", angle::angleIndices::armorDegrees)
        .value("fuse", angle::angleIndices::fuseRadians)
        .value("fuseD", angle::angleIndices::fuseDegrees);

    emscripten::enum_<dispersion::dispersionIndices>("dispersionIndices")
        .value("maxHorizontal", dispersion::dispersionIndices::maxHorizontal)
        .value("standardHorizontal",
               dispersion::dispersionIndices::standardHorizontal)
        .value("halfHorizontal", dispersion::dispersionIndices::halfHorizontal)
        .value("maxVertical", dispersion::dispersionIndices::maxVertical)
        .value("standardVertical",
               dispersion::dispersionIndices::standardVertical)
        .value("halfVertical", dispersion::dispersionIndices::halfVertical)
        .value("maxArea", dispersion::dispersionIndices::maxArea)
        .value("standardArea", dispersion::dispersionIndices::standardArea)
        .value("halfArea", dispersion::dispersionIndices::halfArea);

    emscripten::enum_<post::postPenIndices>("postPenIndices")
        .value("angle", post::postPenIndices::angle)
        .value("distance", post::postPenIndices::distance)
        .value("x", post::postPenIndices::x)
        .value("y", post::postPenIndices::y)
        .value("z", post::postPenIndices::z)
        .value("xwf", post::postPenIndices::xwf);

    emscripten::enum_<calculateType::calcIndices>("calcIndices")
        .value("impact", calculateType::calcIndices::impact)
        .value("angle", calculateType::calcIndices::angle)
        .value("dispersion", calculateType::calcIndices::dispersion)
        .value("post", calculateType::calcIndices::post);
};
}  // namespace wows_shell