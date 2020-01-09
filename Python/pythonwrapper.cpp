//pybind11
#define strdup _strdup //if windows.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <utility>
#include <algorithm>
#include <cstddef>
#include "../shellCPP.hpp"

class shellCombined{
    private:
    shell::shellCalc calc;
    shell::shell s;

    public:
    shellCombined(const double v0, const double caliber, const double krupp, const double mass,
    const double normalization, const double cD, const std::string& name, const double threshold, const double fuseTime){
        s.setValues(v0, caliber, krupp, mass, normalization, cD, name, threshold, fuseTime);
    }

    void calcImpact(){
        calc.calculateImpact(s, false);
    }

    void calcPostPen(const double thickness, const double inclination, std::vector<double> angles){
        calc.calculatePostPen(thickness, inclination, s, angles);
    }

    void printImpact(){
        if(s.completedImpact){
            s.printImpactData();
        }else{
            throw std::runtime_error("Impact data not generated");
        }
    }

    void printPostPen(){
        if(s.completedPostPen){
            s.printPostPenData();
        }else{
            throw std::runtime_error("PostPen data not generated");
        }
    }

    pybind11::array_t<double> getImpact(){
        if(s.completedImpact){
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.getImpactPtr(0, 0),                       /* Pointer to data (nullptr -> ask NumPy to allocate!) */
                sT,                                         /* Size of one item */
                pybind11::format_descriptor<double>::value, /* Buffer format */
                2,                                          /* How many dimensions? */
                std::vector<std::size_t>{ shell::impact::maxColumns, s.impactSize }, /* Number of elements for each dimension */
                std::vector<std::size_t>{ s.impactSizeAligned * sT, sT}                     /* Strides for each dimension */
            ));
            return result;
        }else{
            throw std::runtime_error("Impact data not generated");
        }
    }


    pybind11::array_t<double> getPostPen(){
        if(s.completedPostPen){
            constexpr std::size_t sT = sizeof(double);
            auto result = pybind11::array(pybind11::buffer_info(
                s.getPostPenPtr(0, 0),                      
                sT,                                         
                pybind11::format_descriptor<double>::value, 
                3,                                          
                std::vector<std::size_t>{ shell::post::maxColumns, (int)s.postPenSize / s.impactSize, s.impactSize },  
                std::vector<std::size_t>{ s.postPenSize * sT, s.impactSize * sT, sT}                    
            ));
            return result;
        }else{
            throw std::runtime_error("PostPen data not generated");
        }
    }

};

//Testline: s = shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76.0, .033 )
PYBIND11_MODULE(pythonwrapper, m){
    pybind11::class_<shellCombined>(m, "shell", pybind11::buffer_protocol())
        .def(pybind11::init<double, double, double, double, double, double, std::string& , double, double>())
        .def("calcImpact"  , &shellCombined::calcImpact)
        .def("calcPostPen" , &shellCombined::calcPostPen)
        .def("getImpact"   , &shellCombined::getImpact, pybind11::return_value_policy::reference)
        .def("getPostPen"  , &shellCombined::getPostPen, pybind11::return_value_policy::reference)
        .def("printImpact" , &shellCombined::printImpact)
        .def("printPostPen", &shellCombined::printPostPen);
    
    //Enums
    pybind11::enum_<shell::impact::impactDataIndex>(m, "impactDataIndex", pybind11::arithmetic())
        .value("distance"  , shell::impact::impactDataIndex::distance)
        .value("launchA"   , shell::impact::impactDataIndex::launchA)
        .value("impactAHR" , shell::impact::impactDataIndex::impactAHR)
        .value("impactAHD" , shell::impact::impactDataIndex::impactAHD)
        .value("impactV"   , shell::impact::impactDataIndex::impactV)
        .value("rawPen"    , shell::impact::impactDataIndex::rawPen)
        .value("ePenH"     , shell::impact::impactDataIndex::ePenH)
        .value("ePenHN"    , shell::impact::impactDataIndex::ePenHN)
        .value("impactADD" , shell::impact::impactDataIndex::impactADD)
        .value("ePenD"     , shell::impact::impactDataIndex::ePenD)
        .value("ePenDN"    , shell::impact::impactDataIndex::ePenDN)
        .value("tToTarget" , shell::impact::impactDataIndex::tToTarget)
        .value("tToTargetA", shell::impact::impactDataIndex::tToTargetA)
        .export_values();
    
    pybind11::enum_<shell::post::postPenDataIndex>(m, "postPenDataIndex", pybind11::arithmetic())
        .value("angle"   , shell::post::postPenDataIndex::angle)
        .value("distance", shell::post::postPenDataIndex::distance)
        .value("x"       , shell::post::postPenDataIndex::x)
        .value("y"       , shell::post::postPenDataIndex::y)
        .value("z"       , shell::post::postPenDataIndex::z)
        .value("xwf"     , shell::post::postPenDataIndex::xwf)
        .export_values();
};

