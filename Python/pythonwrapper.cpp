#define strdup _strdup
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <utility>
#include <cstddef>
#include "../shellCPP.hpp"

class shellCombined{
    private:
    shellCalc calc;
    shell s;

    public:
    shellCombined(const double v0, const double caliber, const double krupp, const double mass,
    const double normalization, const double cD, const std::string& name, const double threshold, const double fuseTime){
        s.setValues(v0, caliber, krupp, mass, normalization, cD, name, threshold, fuseTime);
    }

    void calcStandard(){
        calc.calculateStd(s);
        s.printStdData();
    }

    double* stdPtr(){
        return s.stdData.data();
    }

    double* postPenPtr(){
        return s.postPenData.data();
    }

    void calcPostPen(double thickness, std::vector<double> angles){
        s.angles = std::move(angles);
        calc.calculatePostPen(thickness, s);
        s.printPostPen();
    }

    pybind11::array_t<double> getStd(){
        constexpr std::size_t sT = sizeof(double);
        auto result = pybind11::array(pybind11::buffer_info(
            stdPtr(),                                   /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sT,                             /* Size of one item */
            pybind11::format_descriptor<double>::value, /* Buffer format */
            2,                                          /* How many dimensions? */
            std::vector<std::size_t>{ 14, s.sizeAligned },                            /* Number of elements for each dimension */
            std::vector<std::size_t>{ s.sizeAligned * sT, sT}                          /* Strides for each dimension */
        ));
        return result;
    }

    pybind11::array_t<double> getPostPen(){
        constexpr std::size_t sT = sizeof(double);
        auto result = pybind11::array(pybind11::buffer_info(
            postPenPtr(),                       /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sT,                             /* Size of one item */
            pybind11::format_descriptor<double>::value, /* Buffer format */
            3,                                  /* How many dimensions? */
            std::vector<std::size_t>{ s.angles.size(), 6, s.sizeAligned },                          /* Number of elements for each dimension */
            std::vector<std::size_t>{ 6 * s.sizeAligned * sT, s.sizeAligned * sT, sT}                          /* Strides for each dimension */
        ));
        return result;
    }

};

//Testline: s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
PYBIND11_MODULE(pythonwrapper, m){
    pybind11::class_<shellCombined>(m, "shell", pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double, std::string& , double, double>())
        .def("calcStandard", &shellCombined::calcStandard)
        .def("calcPostPen", &shellCombined::calcPostPen)
        .def("getStandard", &shellCombined::getStd)
        .def("getPostPen", &shellCombined::getPostPen);
};
