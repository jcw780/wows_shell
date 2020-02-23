#include "../shellCPP.hpp"
#include <algorithm>
#include <cstddef>
#include <emscripten/bind.h>
#include <utility>

class shellCombined {
private:
    shell::shellCalc calc;
    shell::shell s;

public:
    shellCombined(const double caliber, const double v0, const double cD,
                  const double mass, const double krupp,
                  const double normalization, const double fuseTime,
                  const double threshold, const double ricochet0,
                  const double ricochet1) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, "ship");
    }

    void calcImpact() {
        calc.calculateImpact(s, false, 1); // atomics don't work yet
    }

    int impactSize() { return s.impactSize; }

    int impactSizeAligned() { return s.impactSizeAligned; }

    std::vector<double> impactData() {
        if (s.completedImpact) {
            return s.impactData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getImpactPoint(const int i, const int j) {
        return s.get_impact(i, j);
    }

    void calcPostPen(const double thickness, const double inclination,
                     emscripten::val v, const bool changeDirection,
                     const bool fast) {
        std::vector<double> input =
            std::move(emscripten::vecFromJSArray<double>(v));
        calc.calculatePostPen(thickness, inclination, s, input, changeDirection,
                              fast, 1); // atomics don't work yet
    }

    int postPenSize() { return s.postPenSize; }

    std::vector<double> postPenData() {
        if (s.completedPostPen) {
            return s.postPenData;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    double getPostPenPoint(const int i, const int j, const int k) {
        return s.get_postPen(i, j, k);
    }

    void printImpact() {
        if (s.completedImpact) {
            s.printImpactData();
        } else {
            throw std::runtime_error("Impact data not generated");
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

// Compile option
// emcc --bind -o shellWasm.js shellWasm.cpp --std=c++17 -O3 -s ASSERTIONS=1 -s
// ALLOW_MEMORY_GROWTH=1 emcc --bind -o shellWasm.js shellWasm.cpp --std=c++17
// -O3 -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 -s USE_PTHREADS=1 -s
// WASM_MEM_MAX=100Mb

// Testline: s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
EMSCRIPTEN_BINDINGS(shellWasm) {
    emscripten::class_<shellCombined>("shell")
        .constructor<double, double, double, double, double, double,
                     /*std::string&,*/ double, double>()
        .function("calcImpact", &shellCombined::calcImpact)
        .function("getImpactPoint", &shellCombined::getImpactPoint)
        .function("impactData", &shellCombined::impactData)
        .function("getImpactSize", &shellCombined::impactSize)
        .function("getImpactSizeAligned", &shellCombined::impactSizeAligned)
        .function("calcPostPen", &shellCombined::calcPostPen)
        .function("getPostPenPoint", &shellCombined::getPostPenPoint)
        .function("postPenData", &shellCombined::postPenData)
        .function("getPostPenSize", &shellCombined::postPenSize)
        .function("printImpact", &shellCombined::printImpact)
        .function("printPostPen", &shellCombined::printPostPen);

    emscripten::register_vector<double>("vector<double>");

    // Enums
    emscripten::enum_<shell::impact::impactDataIndex>("impactDataIndex")
        .value("distance", shell::impact::impactDataIndex::distance)
        .value("launchA", shell::impact::impactDataIndex::launchA)
        .value("impactAHR", shell::impact::impactDataIndex::impactAHR)
        .value("impactAHD", shell::impact::impactDataIndex::impactAHD)
        .value("impactV", shell::impact::impactDataIndex::impactV)
        .value("rawPen", shell::impact::impactDataIndex::rawPen)
        .value("ePenH", shell::impact::impactDataIndex::ePenH)
        .value("ePenHN", shell::impact::impactDataIndex::ePenHN)
        .value("impactADD", shell::impact::impactDataIndex::impactADD)
        .value("ePenD", shell::impact::impactDataIndex::ePenD)
        .value("ePenDN", shell::impact::impactDataIndex::ePenDN)
        .value("tToTarget", shell::impact::impactDataIndex::tToTarget)
        .value("tToTargetA", shell::impact::impactDataIndex::tToTargetA);

    emscripten::enum_<shell::post::postPenDataIndex>("postPenDataIndex")
        .value("angle", shell::post::postPenDataIndex::angle)
        .value("distance", shell::post::postPenDataIndex::distance)
        .value("x", shell::post::postPenDataIndex::x)
        .value("y", shell::post::postPenDataIndex::y)
        .value("z", shell::post::postPenDataIndex::z)
        .value("xwf", shell::post::postPenDataIndex::xwf);
};