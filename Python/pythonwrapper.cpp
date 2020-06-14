/*cppimport
<%
cfg['compiler_args'] = ['-std=c++17', '/std:c++17', '-O3', '/Ot',
'-march=native','/arch:AVX2']
setup_pybind11(cfg)
%>
*/

/*
strdup not defined in windows
https://github.com/pybind/pybind11/issues/1212
*/
#ifdef _WIN32
#define strdup _strdup
#endif

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../shellCPP.hpp"
#include <algorithm>
#include <cstddef>
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
                  const double ricochet1, const std::string &name) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, name);
    }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const std::string &name) {
        s.setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                    threshold, ricochet0, ricochet1, name);
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
        calc.calculateImpact<shell::numerical::forwardEuler, false>(s, false);
    }

    void calcImpactAdamsBashforth5() {
        calc.calculateImpact<shell::numerical::adamsBashforth5, false>(s,
                                                                       false);
    }

    void calcImpactRungeKutta2() {
        calc.calculateImpact<shell::numerical::rungeKutta2, false>(s, false);
    }

    void calcImpactRungeKutta4() {
        calc.calculateImpact<shell::numerical::rungeKutta4, false>(s, false);
    }

    void calcImpactRungeKutta4Hybrid() {
        calc.calculateImpact<shell::numerical::rungeKutta4, true>(s, true);
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

    pybind11::array_t<double> getImpact() {
        if (s.completedImpact) {
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_impactPtr(0, 0), /* Pointer to data (nullptr -> ask NumPy
                                         to allocate!) */
                sT,                    /* Size of one item */
                pybind11::format_descriptor<double>::value, /* Buffer format */
                2, /* How many dimensions? */
                std::vector<std::size_t>{
                    shell::impact::maxColumns,
                    s.impactSize}, /* Number of elements for each dimension */
                std::vector<std::size_t>{s.impactSizeAligned * sT, sT}
                /* Strides for each dimension */
                ));
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getAngles() {
        if (s.completedImpact) {
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_anglePtr(0, 0), /* Pointer to data (nullptr -> ask NumPy
                                         to allocate!) */
                sT,                   /* Size of one item */
                pybind11::format_descriptor<double>::value, /* Buffer format */
                2, /* How many dimensions? */
                std::vector<std::size_t>{
                    shell::angle::maxColumns,
                    s.impactSize}, /* Number of elements for each dimension */
                std::vector<std::size_t>{s.impactSizeAligned * sT, sT}
                /* Strides for each dimension */
                ));
            return result;
        } else {
            throw std::runtime_error("Impact data not generated");
        }
    }

    pybind11::array_t<double> getPostPen() {
        if (s.completedPostPen) {
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.get_postPenPtr(0, 0, 0), sT,
                pybind11::format_descriptor<double>::value, 3,
                std::vector<std::size_t>{shell::post::maxColumns,
                                         (int)s.postPenSize / s.impactSize,
                                         s.impactSize},
                std::vector<std::size_t>{s.postPenSize * sT, s.impactSize * sT,
                                         sT}));
            return result;
        } else {
            throw std::runtime_error("PostPen data not generated");
        }
    }
};

PYBIND11_MODULE(pythonwrapper, m) {
    pybind11::class_<shellCombined>(m, "shell", pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double, std::string &>())
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
        .def("calcImpactRungeKutta4Hybrid",
             &shellCombined::calcImpactRungeKutta4Hybrid)
        .def("calcAngles", &shellCombined::calcAngles)
        .def("calcPostPen", &shellCombined::calcPostPen)
        .def("getImpact", &shellCombined::getImpact)
        // pybind11::return_value_policy::reference)
        .def("getAngles", &shellCombined::getAngles)
        .def("getPostPen", &shellCombined::getPostPen)
        // pybind11::return_value_policy::reference)
        .def("printImpact", &shellCombined::printImpact)
        .def("printAngles", &shellCombined::printAngles)
        .def("printPostPen", &shellCombined::printPostPen);

    // Enums
    pybind11::enum_<shell::impact::impactDataIndex>(m, "impactDataIndex",
                                                    pybind11::arithmetic())
        .value("distance", shell::impact::impactDataIndex::distance)
        .value("launchAngle", shell::impact::impactDataIndex::launchAngle)
        .value("impactAngleHorizontalRadians",
               shell::impact::impactDataIndex::impactAngleHorizontalRadians)
        .value("impactAngleHorizontalDegrees",
               shell::impact::impactDataIndex::impactAngleHorizontalDegrees)
        .value("impactVelocity", shell::impact::impactDataIndex::impactVelocity)
        .value("rawPenetration", shell::impact::impactDataIndex::rawPenetration)
        .value("effectivePenetrationHorizontal",
               shell::impact::impactDataIndex::effectivePenetrationHorizontal)
        .value("effectivePenetrationHorizontalNormalized",
               shell::impact::impactDataIndex::
                   effectivePenetrationHorizontalNormalized)
        .value("impactAngleDeckDegrees",
               shell::impact::impactDataIndex::impactAngleDeckDegrees)
        .value("effectivePenetrationDeck",
               shell::impact::impactDataIndex::effectivePenetrationDeck)
        .value(
            "effectivePenetrationDeckNormalized",
            shell::impact::impactDataIndex::effectivePenetrationDeckNormalized)
        .value("timeToTarget", shell::impact::impactDataIndex::timeToTarget)
        .value("timeToTargetAdjusted",
               shell::impact::impactDataIndex::timeToTargetAdjusted)
        .export_values();

    pybind11::enum_<shell::angle::angleDataIndex>(m, "angleDataIndex",
                                                  pybind11::arithmetic())
        .value("distance", shell::angle::angleDataIndex::distance)
        .value("ricochetAngle0Radians",
               shell::angle::angleDataIndex::ricochetAngle0Radians)
        .value("ricochetAngle0Degrees",
               shell::angle::angleDataIndex::ricochetAngle0Degrees)
        .value("ricochetAngle1Radians",
               shell::angle::angleDataIndex::ricochetAngle1Radians)
        .value("ricochetAngle1Degrees",
               shell::angle::angleDataIndex::ricochetAngle1Degrees)
        .value("armorRadians", shell::angle::angleDataIndex::armorRadians)
        .value("armorDegrees", shell::angle::angleDataIndex::armorDegrees)
        .value("fuseRadians", shell::angle::angleDataIndex::fuseRadians)
        .value("fuseDegrees", shell::angle::angleDataIndex::fuseDegrees)
        .export_values();

    pybind11::enum_<shell::post::postPenDataIndex>(m, "postPenDataIndex",
                                                   pybind11::arithmetic())
        .value("angle", shell::post::postPenDataIndex::angle)
        .value("distance", shell::post::postPenDataIndex::distance)
        .value("x", shell::post::postPenDataIndex::x)
        .value("y", shell::post::postPenDataIndex::y)
        .value("z", shell::post::postPenDataIndex::z)
        .value("xwf", shell::post::postPenDataIndex::xwf)
        .export_values();
};
