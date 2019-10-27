from libcpp cimport vector

cdef extern from "shellCPP.cpp":
    cdef cppclass shell:
        shell(double v0, double caliber, double krupp, double mass, double normalization, double cD, std::string name, double threshold, double fuseTime) except +
        void calculateStd()

        
        void calculatePostPen()
        double *exposeStdData()
        double *exposePostPenData()

