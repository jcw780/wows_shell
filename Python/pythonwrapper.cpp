#define strdup _strdup
#include <pybind11/pybind11.h>
#include <utility>
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

    calcStandard(){
        calc.calculateStd(s);
        s.printStdData;
    }

    calcPostPen(double thickness, std::vector<double> angles){
        s.angles = std::move(angles);
        calc.calculatePostPen(thickness, s);
        s.printPostPen();
    }
}


PYBIND11_MODULE(pythonwrapper, m){
    py::class_<shellCombined>(m, "shell")
        .def(py::init<double, double, double, double, double, double, std::string& , double, double>())
        .def("calcStandard", );
        .def("calcPostPen", );
}
