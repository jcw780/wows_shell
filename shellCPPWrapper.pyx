from libcpp cimport vector, string
import numpy as np
cimport numpy as np

cdef extern from "shellCPP.cpp":
    cdef cppclass shell:
        shell(double v0, double caliber, double krupp, double mass, double normalization, double cD, string name, double threshold, double fuseTime) except +
        void calculateStd()

        void setAngles(vector<double> *angles)
        #void setAngles()
        void calculatePostPen(double thickness)
        double *exposeStdData()
        double *exposePostPenData()

        unsigned int size, sizeAligned, postPenSize, angleSize
        

cdef class pyShell:
    cdef shell s
    def __cinit__(self, input):
        if type(input) == dict:
            self.s(input['v0'], input['caliber'], input['krupp'], input['mass'], input['normalization'], input['cD'], input['name'], input['threshold'], input['fuseTime'],)
        else:
            self.s(input[0], input[1], input[2], input[3], input[4], input[5], <string> input[6].encode('utf-8'), input[7], input[8])
        #self.name = input[6]
    
    def calculateStd():
        self.s.calculateStd()
    
    def getStandard():
        return np.asarray(<np.float64[:13*self.s.size]> self.s.exposeStdData())

    def calculatePostPen(thickness):
        self.s.calculatePostPen(thickness)

    def getPostPen():
        return np.asarray(<np.float64[:6*self.s.postPenSize]]> self.s.exposePostPenData())
    
    def setAngles(angles):
        self.s.setAngles(angles)





