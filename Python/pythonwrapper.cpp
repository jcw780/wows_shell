#define strdup _strdup
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <utility>
#include <algorithm>
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
        pybind11::gil_scoped_release release;
        calc.calculateStd(s, false);
        //s.printStdData();
        pybind11::gil_scoped_acquire acquire;
    }

    void calcPostPen(double thickness, std::vector<double> angles){
        //s.angles = std::move(angles);
        std::cout<<"Entered"<<std::endl;
        
        /*s.angles.resize(angles.size());
        std::copy(angles.begin(), angles.end(), s.angles.begin());*/
        //s.angles = angles;
        pybind11::gil_scoped_release release;
        //s.angles = std::move(angles);
        std::cout<<"GIL Released"<<std::endl;
        calc.calculatePostPen(thickness, s, angles);
        std::cout<<"GIL Acquired"<<std::endl;
        //s.printPostPen();
        pybind11::gil_scoped_acquire acquire;
    }

    double* stdPtr(){
        return s.stdData.data();
    }

    double* postPenPtr(){
        return s.postPenData.data();
    }

    pybind11::array_t<double> getStd(){
        constexpr std::size_t sT = sizeof(double);
        //double *temp = new double[s.sizeAligned * 13];
        //std::copy_n(stdPtr(), s.sizeAligned * 13, temp);

        auto result = pybind11::array(pybind11::buffer_info(
            nullptr,                                   /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sT,                             /* Size of one item */
            pybind11::format_descriptor<double>::value, /* Buffer format */
            2,                                          /* How many dimensions? */
            std::vector<std::size_t>{ 13, s.sizeAligned },                            /* Number of elements for each dimension */
            std::vector<std::size_t>{ s.sizeAligned * sT, sT}                          /* Strides for each dimension */
        ));
        std::copy_n(stdPtr(), s.sizeAligned * 13, (double*) result.request().ptr);
        return result;
    }

    pybind11::array_t<double> getPostPen(){
        std::cout<<"Returning"<<std::endl;
        constexpr std::size_t sT = sizeof(double);
        //double *temp = new double[s.postPenSize * 6];
        auto result = pybind11::array(pybind11::buffer_info(
            nullptr,                       /* Pointer to data (nullptr -> ask NumPy to allocate!) */
            sT,                             /* Size of one item */
            pybind11::format_descriptor<double>::value, /* Buffer format */
            2,                                  /* How many dimensions? */
            std::vector<std::size_t>{ 6, s.postPenSize },                          /* Number of elements for each dimension */
            std::vector<std::size_t>{ s.postPenSize * sT, sT}                          /* Strides for each dimension */
        ));
        std::copy_n(postPenPtr(), s.postPenSize * 6, (double*) result.request().ptr);
        std::cout<<"Returning"<<std::endl;
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
