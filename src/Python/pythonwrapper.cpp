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
#ifdef _WIN32
#define strdup _strdup
#endif

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstddef>
#include <utility>

#include "../shellCPP.hpp"

class shellPython {
   public:
    wows_shell::shell s;
    shellPython(const double caliber, const double v0, const double cD,
                const double mass, const double krupp,
                const double normalization, const double fuseTime,
                const double threshold, const double ricochet0,
                const double ricochet1, const double nonAP,
                const std::string &name) {
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
                s.get_impactPtr(0, 0), /* Pointer to data (nullptr -> ask NumPy
                                         to allocate!) */
                sT,                    /* Size of one item */
                pybind11::format_descriptor<double>::value, /* Buffer format */
                2, /* How many dimensions? */
                std::vector<std::size_t>{
                    wows_shell::impact::maxColumns,
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
                    wows_shell::angle::maxColumns,
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
                std::vector<std::size_t>{wows_shell::post::maxColumns,
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

std::string generateShellPythonHash (const shellPython& s){
    return wows_shell::generateHash(s.s);
}

class shellCalcPython : public wows_shell::shellCalc {
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
    
    template <wows_shell::numerical Numerical>
    void calcImpact(shellPython &sp){
        calculateImpact<false, Numerical, false>(sp.s);
    }

    void calcAngles(shellPython &sp, const double thickness,
                    const double inclination) {
        calculateAngles(thickness, inclination, sp.s);
    }

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
    wows_shell::shellCalc calc;
    wows_shell::shell s;

   public:
    shellCombined(const double caliber, const double v0, const double cD,
                  const double mass, const double krupp,
                  const double normalization, const double fuseTime,
                  const double threshold, const double ricochet0,
                  const double ricochet1, const double nonAP,
                  const std::string &name) {
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
        calc.calculateImpact<false, wows_shell::numerical::forwardEuler, false>(s);
    }

    void calcImpactAdamsBashforth5() {
        calc.calculateImpact<false, wows_shell::numerical::adamsBashforth5, false>(
            s);
    }

    void calcImpactRungeKutta2() {
        calc.calculateImpact<false, wows_shell::numerical::rungeKutta2, false>(s);
    }

    void calcImpactRungeKutta4() {
        calc.calculateImpact<false, wows_shell::numerical::rungeKutta4, false>(s);
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
                s.get_impactPtr(0, 0), /* Pointer to data (nullptr -> ask NumPy
                                         to allocate!) */
                sT,                    /* Size of one item */
                pybind11::format_descriptor<double>::value, /* Buffer format */
                2, /* How many dimensions? */
                std::vector<std::size_t>{
                    wows_shell::impact::maxColumns,
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
                    wows_shell::angle::maxColumns,
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
                std::vector<std::size_t>{wows_shell::post::maxColumns,
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
        // pybind11::return_value_policy::reference)
        .def("getAngles", &shellCombined::getAngles)
        .def("getPostPen", &shellCombined::getPostPen)
        // pybind11::return_value_policy::reference)
        .def("printImpact", &shellCombined::printImpact)
        .def("printAngles", &shellCombined::printAngles)
        .def("printPostPen", &shellCombined::printPostPen);
    
    m.def("generateHash", &wows_shell::generateShellParamHash);
    m.def("generateShellHash", &generateShellPythonHash);

    pybind11::class_<shellPython>(m, "shell", pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double,
                            double, double, double, double, double,
                            std::string &>())
        .def("setValues", &shellPython::setValues)
        .def("interpolateDistanceImpact",
             &shellPython::interpolateDistanceImpact)
        .def("getImpact", &shellPython::getImpact)
        // pybind11::return_value_policy::reference)
        .def("getAngles", &shellPython::getAngles)
        .def("getPostPen", &shellPython::getPostPen)
        // pybind11::return_value_policy::reference)
        .def("printImpact", &shellPython::printImpact)
        .def("printAngles", &shellPython::printAngles)
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
        .def("calcImpactForwardEuler", &shellCalcPython::calcImpact<wows_shell::numerical::forwardEuler>)
        .def("calcImpactAdamsBashforth5",
             &shellCalcPython::calcImpact<wows_shell::numerical::adamsBashforth5>)
        .def("calcImpactRungeKutta2", &shellCalcPython::calcImpact<wows_shell::numerical::rungeKutta2>)
        .def("calcImpactRungeKutta4", &shellCalcPython::calcImpact<wows_shell::numerical::rungeKutta4>)
        .def("calcAngles", &shellCalcPython::calcAngles)
        .def("calcPostPen", &shellCalcPython::calcPostPen);
    // Enums
    pybind11::enum_<wows_shell::impact::impactIndices>(m, "impactIndices",
                                                  pybind11::arithmetic())
        .value("distance", wows_shell::impact::impactIndices::distance)
        .value("launchAngle", wows_shell::impact::impactIndices::launchAngle)
        .value("impactAngleHorizontalRadians",
               wows_shell::impact::impactIndices::impactAngleHorizontalRadians)
        .value("impactAngleHorizontalDegrees",
               wows_shell::impact::impactIndices::impactAngleHorizontalDegrees)
        .value("impactVelocity", wows_shell::impact::impactIndices::impactVelocity)
        .value("rawPenetration", wows_shell::impact::impactIndices::rawPenetration)
        .value("effectivePenetrationHorizontal",
               wows_shell::impact::impactIndices::effectivePenetrationHorizontal)
        .value("effectivePenetrationHorizontalNormalized",
               wows_shell::impact::impactIndices::
                   effectivePenetrationHorizontalNormalized)
        .value("impactAngleDeckDegrees",
               wows_shell::impact::impactIndices::impactAngleDeckDegrees)
        .value("effectivePenetrationDeck",
               wows_shell::impact::impactIndices::effectivePenetrationDeck)
        .value("effectivePenetrationDeckNormalized",
               wows_shell::impact::impactIndices::effectivePenetrationDeckNormalized)
        .value("timeToTarget", wows_shell::impact::impactIndices::timeToTarget)
        .value("timeToTargetAdjusted",
               wows_shell::impact::impactIndices::timeToTargetAdjusted)
        .export_values();

    pybind11::enum_<wows_shell::angle::angleIndices>(m, "angleIndices",
                                                pybind11::arithmetic())
        .value("distance", wows_shell::angle::angleIndices::distance)
        .value("ricochetAngle0Radians",
               wows_shell::angle::angleIndices::ricochetAngle0Radians)
        .value("ricochetAngle0Degrees",
               wows_shell::angle::angleIndices::ricochetAngle0Degrees)
        .value("ricochetAngle1Radians",
               wows_shell::angle::angleIndices::ricochetAngle1Radians)
        .value("ricochetAngle1Degrees",
               wows_shell::angle::angleIndices::ricochetAngle1Degrees)
        .value("armorRadians", wows_shell::angle::angleIndices::armorRadians)
        .value("armorDegrees", wows_shell::angle::angleIndices::armorDegrees)
        .value("fuseRadians", wows_shell::angle::angleIndices::fuseRadians)
        .value("fuseDegrees", wows_shell::angle::angleIndices::fuseDegrees)
        .export_values();

    pybind11::enum_<wows_shell::post::postPenIndices>(m, "postPenIndices",
                                                 pybind11::arithmetic())
        .value("angle", wows_shell::post::postPenIndices::angle)
        .value("distance", wows_shell::post::postPenIndices::distance)
        .value("x", wows_shell::post::postPenIndices::x)
        .value("y", wows_shell::post::postPenIndices::y)
        .value("z", wows_shell::post::postPenIndices::z)
        .value("xwf", wows_shell::post::postPenIndices::xwf)
        .export_values();
};
