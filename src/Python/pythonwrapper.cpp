/*cppimport
<%
cfg['compiler_args'] = ['-std=c++17', '/std:c++17', '-Ofast', '/Ot',
'-march=native','/arch:AVX2', '/D_USE_MATH_DEFINES', '/GL']
setup_pybind11(cfg)
%>
*/

/*
strdup not defined in windows
https://github.com/pybind/pybind11/issues/1212
*/
#include <cstdlib>
#include <iterator>
#ifdef _WIN32
#define strdup _strdup
#endif

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <utility>

#include "../shellCPP.hpp"

namespace wows_shell {
class shellPython {
   public:
    shell s;
    shellPython(const double caliber, const double v0, const double cD,
                const double mass, const double krupp,
                const double normalization, const double fuseTime,
                const double threshold, const double ricochet0,
                const double ricochet1, const double nonAP,
                const std::string &name) {
        PyErr_WarnEx(PyExc_DeprecationWarning,
                     "shell(args...) is deprecated, use "
                     "shell(shellParams(args...), name) instead.",
                     1);
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, nonAP, name);
    }
    shellPython(const shellParams &sp, const std::string &name) {
        s.setValues(sp, name);
    }
    shellPython(const shellParams &sp, const dispersionParams &dp,
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

    void printDispersion() {
        if (s.completedDispersion) {
            s.printDispersionData();
        } else {
            throw std::runtime_error("Dispersion data not generated");
        }
    }

    void printPostPen() {
        if (s.completedPostPen) {
            s.printPostPenData();
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }

    double interpolateDistanceImpact(double distance, unsigned int impact) {
        return s.interpolateDistanceImpact(distance, impact);
    }

    pybind11::array_t<double> getImpact(bool owned = true) {
        if (s.completedImpact) {
            constexpr std::size_t sT = sizeof(double);
            std::array<size_t, 2> shape = {impact::maxColumns, s.impactSize},
                                  stride = {s.impactSizeAligned * sT, sT};
            double *tgt = s.get_impactPtr(0, 0);
            auto result =
                owned ? pybind11::array_t<double>(pybind11::buffer_info(
                            tgt, sT, pybind11::format_descriptor<double>::value,
                            2, shape, stride))
                      : pybind11::array_t<double>(shape, stride, tgt);
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getAngles(bool owned = true) {
        if (s.completedAngles) {
            constexpr std::size_t sT = sizeof(double);
            std::array<size_t, 2> shape = {angle::maxColumns, s.impactSize},
                                  stride = {s.impactSizeAligned * sT, sT};
            double *tgt = s.get_anglePtr(0, 0);
            auto result =
                owned ? pybind11::array_t<double>(pybind11::buffer_info(
                            tgt, sT, pybind11::format_descriptor<double>::value,
                            2, shape, stride))
                      : pybind11::array_t<double>(shape, stride, tgt);
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getDispersion(bool owned = true) {
        if (s.completedDispersion) {
            constexpr std::size_t sT = sizeof(double);
            std::array<size_t, 2> shape = {dispersion::maxColumns,
                                           s.impactSize},
                                  stride = {s.impactSizeAligned * sT, sT};
            double *tgt = s.get_dispersionPtr(0, 0);
            auto result =
                owned ? pybind11::array_t<double>(pybind11::buffer_info(
                            tgt, sT, pybind11::format_descriptor<double>::value,
                            2, shape, stride))
                      : pybind11::array_t<double>(shape, stride, tgt);
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getPostPen(bool owned = true) {
        if (s.completedPostPen) {
            constexpr std::size_t sT = sizeof(double);
            std::size_t numAngles = s.postPenSize / s.impactSize;
            std::array<size_t, 3> shape = {post::maxColumns, numAngles,
                                           s.impactSize},
                                  stride = {s.postPenSize * sT,
                                            s.impactSize * sT, sT};
            double *tgt = s.get_postPenPtr(0, 0, 0);
            auto result =
                owned ? pybind11::array_t<double>(pybind11::buffer_info(
                            tgt, sT, pybind11::format_descriptor<double>::value,
                            3, shape, stride))
                      : pybind11::array_t<double>(shape, stride, tgt);
            return result;
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};

std::string generateShellPythonHash(const shellPython &s) {
    return generateHash(s.s);
}

class shellCalcPython : public shellCalc {
   public:
    shellCalcPython() = default;

    /*void setMax(const double max) { calc.set_max(max); }
    void setMin(const double min) { calc.set_min(min); }
    void setPrecision(const double precision) { calc.set_precision(precision); }
    void setX0(const double x0) { calc.set_x0(x0); }
    void setY0(const double y0) { calc.set_y0(y0); }
    void setDtMin(const double dt) { calc.set_dt_min(dt); }
    void setXf0(const double xf0) { calc.set_xf0(xf0); }
    void setYf0(const double yf0) { calc.set_yf0(yf0); }
    void setDtf(const double dtf) { calc.set_dtf(dtf); }*/

    template <numerical Numerical>
    void calcImpact(shellPython &sp) {
        calculateImpact<false, Numerical, false>(sp.s);
    }

    void calcAngles(shellPython &sp, const double thickness,
                    const double inclination) {
        calculateAngles(thickness, inclination, sp.s);
    }

    void calcDispersion(shellPython &sp) { calculateDispersion(sp.s); }

    void calcPostPen(shellPython &sp, const double thickness,
                     const double inclination, std::vector<double> angles,
                     const bool changeDirection, const bool fast) {
        calculatePostPen(thickness, inclination, sp.s, angles, changeDirection,
                         fast);
    }
};

// DEPRECATED
class shellCombined {
   private:
    shellCalc calc;
    shell s;

   public:
    shellCombined(const double caliber, const double v0, const double cD,
                  const double mass, const double krupp,
                  const double normalization, const double fuseTime,
                  const double threshold, const double ricochet0,
                  const double ricochet1, const double nonAP,
                  const std::string &name) {
        PyErr_WarnEx(
            PyExc_DeprecationWarning,
            "shellCombined is deprecated, use shell and shellCalc instead.", 1);
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, nonAP, name);
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

    void setMax(const double max) { calc.set_max(max); }
    void setMin(const double min) { calc.set_min(min); }
    void setPrecision(const double precision) { calc.set_precision(precision); }
    void setX0(const double x0) { calc.set_x0(x0); }
    void setY0(const double y0) { calc.set_y0(y0); }
    void setDtMin(const double dt) { calc.set_dt_min(dt); }
    void setXf0(const double xf0) { calc.set_xf0(xf0); }
    void setYf0(const double yf0) { calc.set_yf0(yf0); }
    void setDtf(const double dtf) { calc.set_dtf(dtf); }

    void calcImpactForwardEuler() {
        calc.calculateImpact<false, numerical::forwardEuler, false>(s);
    }

    void calcImpactAdamsBashforth5() {
        calc.calculateImpact<false, numerical::adamsBashforth5, false>(s);
    }

    void calcImpactRungeKutta2() {
        calc.calculateImpact<false, numerical::rungeKutta2, false>(s);
    }

    void calcImpactRungeKutta4() {
        calc.calculateImpact<false, numerical::rungeKutta4, false>(s);
    }

    void calcAngles(const double thickness, const double inclination) {
        calc.calculateAngles(thickness, inclination, s);
    }

    void calcPostPen(const double thickness, const double inclination,
                     std::vector<double> angles, const bool changeDirection,
                     const bool fast) {
        calc.calculatePostPen(thickness, inclination, s, angles,
                              changeDirection, fast);
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

    double interpolateDistanceImpact(double distance, unsigned int impact) {
        return s.interpolateDistanceImpact(distance, impact);
    }

    pybind11::array_t<double> getImpact() {
        if (s.completedImpact) {
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_impactPtr(0, 0), sT,
                pybind11::format_descriptor<double>::value, 2,
                {impact::maxColumns, s.impactSize},
                {s.impactSizeAligned * sT, sT}));
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getAngles() {
        if (s.completedAngles) {
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_anglePtr(0, 0), sT,
                pybind11::format_descriptor<double>::value, 2,
                {angle::maxColumns, s.impactSize},
                {s.impactSizeAligned * sT, sT}));
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getPostPen() {
        if (s.completedPostPen) {
            constexpr std::size_t sT = sizeof(double);
            std::size_t numAngles = s.postPenSize / s.impactSize;
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_postPenPtr(0, 0, 0), sT,
                pybind11::format_descriptor<double>::value, 3,
                {post::maxColumns, numAngles, s.impactSize},
                {s.postPenSize * sT, s.impactSize * sT, sT}));
            return result;
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};

template <typename Input, typename Keys, typename Output, typename KeyGenerator>
void extractDictToArray(Input &input, Keys &keys, Output &output,
                        KeyGenerator keyGenerator) {
    using V = typename Output::value_type;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        auto adjKey = keyGenerator(keys[i]);
        if (input.contains(adjKey)) {
            output[i] = input[adjKey].template cast<V>();
            std::cout << "Key: " << keys[i] << " Value: " << output[i]
                      << " Added\n";
        } else {
            std::cout << "Key: " << keys[i] << "Not Found\n";
        }
    }
}

template <typename KV>
std::unique_ptr<shellParams> makeShellParamsFromKV(KV &input) {
    constexpr std::size_t structSize = 11;
    constexpr std::array<char const *, structSize> doubleKeys = {
        "caliber",       "v0",       "cD",        "mass",      "krupp",
        "normalization", "fuseTime", "threshold", "ricochet0", "ricochet1",
        "nonAP"};
    std::array<double, structSize> doubleValues{};
    extractDictToArray(input, doubleKeys, doubleValues,
                       [](const char *in) { return pybind11::str(in); });

    return std::make_unique<shellParams>(
        doubleValues[0], doubleValues[1], doubleValues[2], doubleValues[3],
        doubleValues[4], doubleValues[5], doubleValues[6], doubleValues[7],
        doubleValues[8], doubleValues[9], doubleValues[10]);
}

template <typename KV>
std::unique_ptr<dispersionParams> makeDispersionParamsKV(KV &input) {
    constexpr std::size_t structSize = 10;
    constexpr std::array<char const *, structSize> doubleKeys = {
        "idealRadius", "minRadius",  "idealDistance", "taperDistance",
        "delim",       "zeroRadius", "delimRadius",   "maxRadius",
        "maxDistance", "sigma"};
    std::array<double, structSize> doubleValues{};
    extractDictToArray(input, doubleKeys, doubleValues,
                       [](const char *in) { return pybind11::str(in); });

    return std::make_unique<dispersionParams>(
        doubleValues[0], doubleValues[1], doubleValues[2], doubleValues[3],
        doubleValues[4], doubleValues[5], doubleValues[6], doubleValues[7],
        doubleValues[8], doubleValues[9]);
}

PYBIND11_MODULE(pythonwrapper, m) {
    pybind11::class_<shellParams>(m, "shellParams")
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double, double>())
        .def(pybind11::init(
            [](pybind11::dict &i) { return makeShellParamsFromKV(i); }))
        .def(pybind11::init(
            [](pybind11::kwargs &i) { return makeShellParamsFromKV(i); }))
        .def(pybind11::init())
        .def_readwrite("caliber", &shellParams::caliber)
        .def_readwrite("v0", &shellParams::v0)
        .def_readwrite("cD", &shellParams::cD)
        .def_readwrite("mass", &shellParams::mass)
        .def_readwrite("krupp", &shellParams::krupp)
        .def_readwrite("normalization", &shellParams::normalization)
        .def_readwrite("fuseTime", &shellParams::fuseTime)
        .def_readwrite("threshold", &shellParams::threshold)
        .def_readwrite("ricochet0", &shellParams::ricochet0)
        .def_readwrite("ricochet1", &shellParams::ricochet1)
        .def_readwrite("nonAP", &shellParams::nonAP);

    pybind11::class_<dispersionParams>(m, "dispersionParams")
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double>())
        .def(pybind11::init(
            [](pybind11::dict &i) { return makeDispersionParamsKV(i); }))
        .def(pybind11::init(
            [](pybind11::kwargs &i) { return makeDispersionParamsKV(i); }))
        .def(pybind11::init())
        .def_readwrite("idealRadius", &dispersionParams::idealRadius)
        .def_readwrite("minRadius", &dispersionParams::minRadius)
        .def_readwrite("idealDistance", &dispersionParams::idealDistance)
        .def_readwrite("taperDistance", &dispersionParams::taperDistance)
        .def_readwrite("delim", &dispersionParams::delim)
        .def_readwrite("zeroRadius", &dispersionParams::zeroRadius)
        .def_readwrite("delimRadius", &dispersionParams::delimRadius)
        .def_readwrite("maxRadius", &dispersionParams::maxRadius)
        .def_readwrite("maxDistance", &dispersionParams::maxDistance)
        .def_readwrite("sigma", &dispersionParams::sigma);

    pybind11::class_<shellCombined>(m, "shellCombined",
                                    pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double, double,
                            std::string &>())
        .def("setValues", &shellCombined::setValues)
        .def("setMax", &shellCombined::setMax)
        .def("setMin", &shellCombined::setMin)
        .def("setPrecision", &shellCombined::setPrecision)
        .def("setX0", &shellCombined::setX0)
        .def("setY0", &shellCombined::setY0)
        .def("setDtMin", &shellCombined::setDtMin)
        .def("setXf0", &shellCombined::setXf0)
        .def("setYf0", &shellCombined::setYf0)
        .def("setDtf", &shellCombined::setDtf)

        .def("calcImpactForwardEuler", &shellCombined::calcImpactForwardEuler)
        .def("calcImpactAdamsBashforth5",
             &shellCombined::calcImpactAdamsBashforth5)
        .def("calcImpactRungeKutta2", &shellCombined::calcImpactRungeKutta2)
        .def("calcImpactRungeKutta4", &shellCombined::calcImpactRungeKutta4)

        .def("calcAngles", &shellCombined::calcAngles)
        .def("calcPostPen", &shellCombined::calcPostPen)
        .def("interpolateDistanceImpact",
             &shellCombined::interpolateDistanceImpact)
        .def("getImpact", &shellCombined::getImpact)
        .def("getAngles", &shellCombined::getAngles)
        .def("getPostPen", &shellCombined::getPostPen)
        .def("printImpact", &shellCombined::printImpact)
        .def("printAngles", &shellCombined::printAngles)
        .def("printPostPen", &shellCombined::printPostPen);

    m.def("generateHash", &generateShellParamHash);
    m.def("generateShellHash", &generateShellPythonHash);

    pybind11::class_<shellPython>(m, "shell", pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double, double,
                            std::string &>())
        .def(pybind11::init<shellParams &, std::string &>())
        .def(pybind11::init<shellParams &, dispersionParams &, std::string &>())
        .def(
            "setValues",
            static_cast<void (shellPython::*)(
                const double, const double, const double, const double,
                const double, const double, const double, const double,
                const double, const double, const double, const std::string &)>(
                &shellPython::setValues))
        .def("setValues", static_cast<void (shellPython::*)(
                              const shellParams &, const std::string &)>(
                              &shellPython::setValues))
        .def("setValues", static_cast<void (shellPython::*)(
                              const shellParams &, const dispersionParams &,
                              const std::string &)>(&shellPython::setValues))
        .def("interpolateDistanceImpact",
             &shellPython::interpolateDistanceImpact)
        .def("getImpact", &shellPython::getImpact,
             pybind11::arg("owned") = true)
        .def("getAngles", &shellPython::getAngles,
             pybind11::arg("owned") = true)
        .def("getDispersion", &shellPython::getDispersion,
             pybind11::arg("owned") = true)
        .def("getPostPen", &shellPython::getPostPen,
             pybind11::arg("owned") = true)
        .def("printImpact", &shellPython::printImpact)
        .def("printAngles", &shellPython::printAngles)
        .def("printDispersion", &shellPython::printDispersion)
        .def("printPostPen", &shellPython::printPostPen);

    pybind11::class_<shellCalcPython>(m, "shellCalc",
                                      pybind11::buffer_protocol())
        .def(pybind11::init())
        .def("setMax", &shellCalcPython::set_max)
        .def("setMin", &shellCalcPython::set_min)
        .def("setPrecision", &shellCalcPython::set_precision)
        .def("setX0", &shellCalcPython::set_x0)
        .def("setY0", &shellCalcPython::set_y0)
        .def("setDtMin", &shellCalcPython::set_dt_min)
        .def("setXf0", &shellCalcPython::set_xf0)
        .def("setYf0", &shellCalcPython::set_yf0)
        .def("setDtf", &shellCalcPython::set_dtf)
        .def("calcImpactForwardEuler",
             &shellCalcPython::calcImpact<numerical::forwardEuler>)
        .def("calcImpactAdamsBashforth5",
             &shellCalcPython::calcImpact<numerical::adamsBashforth5>)
        .def("calcImpactRungeKutta2",
             &shellCalcPython::calcImpact<numerical::rungeKutta2>)
        .def("calcImpactRungeKutta4",
             &shellCalcPython::calcImpact<numerical::rungeKutta4>)
        .def("calcAngles", &shellCalcPython::calcAngles)
        .def("calcDispersion", &shellCalcPython::calcDispersion)
        .def("calcPostPen", &shellCalcPython::calcPostPen);
    // Enums
    pybind11::enum_<impact::impactIndices>(m, "impactIndices",
                                           pybind11::arithmetic())
        .value("distance", impact::impactIndices::distance)
        .value("launchAngle", impact::impactIndices::launchAngle)
        .value("impactAngleHorizontalRadians",
               impact::impactIndices::impactAngleHorizontalRadians)
        .value("impactAngleHorizontalDegrees",
               impact::impactIndices::impactAngleHorizontalDegrees)
        .value("impactVelocity", impact::impactIndices::impactVelocity)
        .value("rawPenetration", impact::impactIndices::rawPenetration)
        .value("effectivePenetrationHorizontal",
               impact::impactIndices::effectivePenetrationHorizontal)
        .value("effectivePenetrationHorizontalNormalized",
               impact::impactIndices::effectivePenetrationHorizontalNormalized)
        .value("impactAngleDeckDegrees",
               impact::impactIndices::impactAngleDeckDegrees)
        .value("effectivePenetrationDeck",
               impact::impactIndices::effectivePenetrationDeck)
        .value("effectivePenetrationDeckNormalized",
               impact::impactIndices::effectivePenetrationDeckNormalized)
        .value("timeToTarget", impact::impactIndices::timeToTarget)
        .value("timeToTargetAdjusted",
               impact::impactIndices::timeToTargetAdjusted)
        .export_values();

    pybind11::enum_<angle::angleIndices>(m, "angleIndices",
                                         pybind11::arithmetic())
        .value("distance", angle::angleIndices::distance)
        .value("ricochetAngle0Radians",
               angle::angleIndices::ricochetAngle0Radians)
        .value("ricochetAngle0Degrees",
               angle::angleIndices::ricochetAngle0Degrees)
        .value("ricochetAngle1Radians",
               angle::angleIndices::ricochetAngle1Radians)
        .value("ricochetAngle1Degrees",
               angle::angleIndices::ricochetAngle1Degrees)
        .value("armorRadians", angle::angleIndices::armorRadians)
        .value("armorDegrees", angle::angleIndices::armorDegrees)
        .value("fuseRadians", angle::angleIndices::fuseRadians)
        .value("fuseDegrees", angle::angleIndices::fuseDegrees)
        .export_values();

    pybind11::enum_<post::postPenIndices>(m, "postPenIndices",
                                          pybind11::arithmetic())
        .value("angle", post::postPenIndices::angle)
        .value("distance", post::postPenIndices::distance)
        .value("x", post::postPenIndices::x)
        .value("y", post::postPenIndices::y)
        .value("z", post::postPenIndices::z)
        .value("xwf", post::postPenIndices::xwf)
        .export_values();
};
}  // namespace wows_shell