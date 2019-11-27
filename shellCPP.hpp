#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
#include <omp.h>
//#include <immintrin.h>
//#include <x86intrin.h>
//#include "constant.h"

#define USE_SIMD

#ifdef _SINGLE_PRECISION
typedef float fPType;
#else
typedef double fPType;
#endif


//WARNING: CURRENTLY AVX2+ ONLY
#ifdef USE_SIMD
#define ENABLE_AVX2
#include "sleef-2.80/simd/sleefsimddp.c"
const unsigned int vSize = 32 / sizeof(fPType);

/*
#ifdef _SINGLE_PRECISION
typedef __m256 fVType
#else
typedef __m256d fVType
#endif
*/

#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <cstring>
#include <algorithm>

/*
double operator"" _kg (long double input){return input;}
double operator"" _lbs(long double input){return input * 0.453592;}

double operator"" _mps(long double input){return input;}
double operator"" _fps(long double input){return input * 0.3048;}
*/

typedef struct{
    double v0;
    double caliber;
    double krupp; 
    double mass;
    double cD; 
    double normalization;
    double threshold; 
    double fuseTime;
    //std::string name;

}shellParams;

typedef struct{
    double v0;
    double k;
    double cw_2;
    double pPPC;
    double normalizationR;
    double threshold;
}procCharacteristic; 

class shell{
    private:
    //double v0;
    double caliber;
    double krupp; 
    double mass;
    double cD; 
    double normalization;
    //double threshold; 
    double fuseTime;
    std::string name;

    //double k, cw_2, pPPC, normalizationR;
    procCharacteristic p;

    void preProcess(){
        p.k = 0.5 * cD * pow((caliber/2),2) * M_PI / mass;
        p.cw_2 = 100+1000/3*caliber;
        p.pPPC = 0.5561613 * krupp/2400 * pow(mass,0.55) / pow((caliber*1000),0.65);
        p.normalizationR = normalization / 180 * M_PI;
    }

    public:
    unsigned int size, postPenSize;
    unsigned int sizeAligned, postPenSizeAligned;
    bool completed = false;
    //procCharacteristic p;
    /*trajectories output
    [0           ]trajx 0        [1           ]trajy 1
    ...
    [size * 2 - 2]trajx size - 1 [size * 2 - 1]trajy size - 1
    */
    std::vector<std::vector<double>> trajectories;
    /*standard output - Indicies refer to multiples of size [Num Colums of 12
    [0 : 1)-distance,          [1 : 2)-launch angle,      [2 : 3)-impact angle - R
    [3 : 4)-impact angle - D,  [4 : 5)-impact velocity,   [5 : 6)-raw pen,           
    [6 : 7)-effectivepenH,     [7 : 8)-effectivepenH w/N, [8 : 9)-IA_d,              
    [9 :10)-effectivepenD,     [10:11)-effectivepenD w/N, [11:12)-ttt,          
    [12:13)-ttta
    */
    std::vector<double> stdData;
    std::vector<double> angles; 
    std::vector<double> postPenData;
    const unsigned int maxColumns = 13;
    enum stdDataIndex{
        distance   , launchA, impactAHR, 
        impactAHD  , impactV, rawPen   , 
        ePenH      , ePenHN , impactADD, 
        ePenD      , ePenDN , tToTarget, 
        tToTargetA};

    shell(double v0, double caliber, double krupp, double mass,
    double normalization, double cD, std::string name, double threshold, double fuseTime){
        this->fuseTime = fuseTime;

        //this->v0 = v0;
        p.v0 = v0;
        this->caliber = caliber;
        this->krupp = krupp;
        this->mass = mass;
        this->normalization = normalization;
        this->cD = cD;
        this->name = name;

        if(threshold){
            p.threshold = threshold;
        }else{
            p.threshold = caliber / 6;
        }
        preProcess();
    }
    /*
    shell(shellParams sp, string name){

    }*/

    procCharacteristic retPC(){
        return p;
    }

    shellParams returnShipParams(){
        shellParams ret;
        ret.caliber = caliber;
        ret.cD = cD;
        ret.fuseTime = fuseTime;
        ret.krupp = krupp;
        ret.mass = mass;
        //ret.name = name;
        ret.normalization = normalization;
        ret.threshold = p.threshold;
        ret.v0 = p.v0;
        return ret;
    }
    
    void printPostPen(){
        for(unsigned int i=0; i<postPenSize; i++){
            for(int j=0; j<6; j++){
                //printf("%f ", postPenData[i + j * postPenSize]);
                std::cout<< std::fixed<< std::setprecision(4) << postPenData[i + j * postPenSize] << " "; 
            }
            //printf("\n");
            std::cout<<"\n";
        }
        //printf("Completed Print\n");
        std::cout<<"Completed Print\n";
    }

    void printStdData(){
        for(unsigned int i=0; i<size; i++){
            for(unsigned int j=0; j<maxColumns; j++){
                std::cout<<std::fixed<<std::setprecision(4)<<stdData[i + sizeAligned * j]<< " ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"Completed"<<std::endl;
    }
    void printTrajectory(unsigned int target){
        if(target >= size){
            std::cout<<"Target Not Within Range of: "<< size<< std::endl;
        }else{
            printf("Index:[%d] X Y\n", target);
            for(std::vector<double>::size_type i = 0; i < trajectories[target*2].size(); i++) {
                std::cout<<trajectories[target*2][i]<<" "<<trajectories[target*2+1][i]<<std::endl;
            }
        }
    }

};

class shellCalc{
    private:
    //double C = 0.5561613;
    double g = 9.81;
    double t0 = 288;
    double L = 0.0065;
    double p0 = 101325;
    double R = 8.31447;
    double M = 0.0289644;
    double cw_1 = 1;

    double max = 25;
    double min = 0;
    double precision = .1;
    double x0 = 0, y0 = 0;
    double dt = .01;

    //double k, cw_2, pPPC, normalizationR;

    std::vector<double> oneVector;
    std::vector<double> temp;

    unsigned int alignmentRequired = vSize * 8;
    /*standardOut1 - indicies refer to multiples of size
    [0:1) v_x [1:2) v_y
    */

    void singleTraj(unsigned int i, shell& s, bool addTraj){
        //std::cout<<"Running0 "<< i<<std::endl;
        const procCharacteristic pC = s.retPC();
        double k = pC.k;
        double cw_2 = pC.cw_2;
        double pPPC = pC.pPPC;
        double normalizationR = pC.normalizationR;
        __m256d angleSIMD, angleRSIMD, temp, v0SIMD = _mm256_set1_pd(pC.v0);
        __m256d vx, vy, tSIMD;

        #ifdef __clang__
        temp = {0.0, 1.0, 2.0, 3.0};
        #else
        temp[0] = 0.0;
        temp[1] = 1.0;
        temp[2] = 2.0;
        temp[3] = 3.0;
        #endif

        angleSIMD = _mm256_fmadd_pd(
            _mm256_add_pd(_mm256_set1_pd(i), temp), 
            _mm256_set1_pd(precision),
            _mm256_set1_pd(min));
        angleRSIMD = _mm256_mul_pd(angleSIMD, _mm256_set1_pd((M_PI / 180)));
        _mm256_storeu_pd(s.stdData.data() + i + s.launchA * s.sizeAligned, angleSIMD);

        vx = _mm256_mul_pd(v0SIMD, xcos(angleRSIMD));
        vy = _mm256_mul_pd(v0SIMD, xsin(angleRSIMD));

        double T, p, rho, t, x, y, v_x, v_y;
        unsigned int counter;
        #define __TrajBuffer__ 128
        double xT[__TrajBuffer__], yT[__TrajBuffer__];

        for(unsigned int j = 0; (j+i<s.size) && (j < vSize); j++){

            if(addTraj){
                s.trajectories[2*(i+j)  ].reserve(__TrajBuffer__);
                s.trajectories[2*(i+j)+1].reserve(__TrajBuffer__);
            }
            //double T, p, rho, t, x, y, v_x, v_y;
            v_x = vx[j];
            v_y = vy[j];
            //std::cout<<vx[j]<<" "<<vy[j]<<std::endl;
            //std::cout<<k<<" "<<v_x<<" "<<v_y<<std::endl;
            x = x0;
            y = y0;
            t = 0;
            s.trajectories[2*(i+j)  ].push_back(x);
            s.trajectories[2*(i+j)+1].push_back(y);
            //printf("%f ",dt);

            while(y >= 0){
                for(counter = 0; counter < __TrajBuffer__; counter++){
                    x += dt*v_x;
                    y += dt*v_y;
                    //printf("%f ", x);
                    T = t0 - L*y;
                    p = p0*pow((1-L*y/t0),(g*M/(R*L)));
                    rho = p*M/(R*T);

                    v_x -= dt*k*rho*(cw_1*v_x*v_x+cw_2*v_x);
                    v_y -= dt*(g - k*rho*(cw_1*v_y*v_y+cw_2*fabs(v_y))*signum(v_y));
                    t += dt;
                    xT[counter] = x;
                    yT[counter] = y;
                    //if(i+j == 248 || i+j == 247){
                    //    std::cout<<i+k<<" "<<T<<" "<<p<<" "<<x<<" "<<y<<" "<<v_x<<" "<<v_y<<std::endl; 
                    //}
                    if(y < 0){
                        break;
                    }
                }
                if(addTraj){
                    s.trajectories[2*(i+j)  ].insert(s.trajectories[2*(i+j)  ].end(), xT, &xT[counter]);
                    s.trajectories[2*(i+j)+1].insert(s.trajectories[2*(i+j)+1].end(), yT, &yT[counter]);

                }
            }
            s.stdData[i+j+s.sizeAligned*s.distance] = x;
            vx[j] = v_x;
            vy[j] = v_y;
            tSIMD[j] = t;
            //std::cout<<k<<" "<<v_x<<std::endl;
        }

        //Calculate [2]IA , [7]IA_D
        __m256d iVSIMD, rPSIMD;  
        angleRSIMD = xatan(_mm256_div_pd(vy,vx));

        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.impactAHR, angleRSIMD);
        _mm256_storeu_pd(
            s.stdData.data()+i+s.sizeAligned*s.impactAHD, 
            _mm256_mul_pd(angleRSIMD, _mm256_set1_pd(180 / M_PI))
        );

        //iADRSIMD = float64_vec(M_PI / 2) + iARSIMD;
        angleSIMD = _mm256_add_pd(
            _mm256_set1_pd(M_PI / 2), angleRSIMD
        );

        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.impactADD,
            _mm256_mul_pd(angleSIMD, _mm256_set1_pd(180 / M_PI))
        );

        //Calculate [3]iV,  [4]rP

        iVSIMD = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(vx, vx), _mm256_mul_pd(vy, vy)));
        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.impactV, iVSIMD);

        rPSIMD = _mm256_mul_pd(xpow(iVSIMD, _mm256_set1_pd(1.1)), _mm256_set1_pd(pPPC));
        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.rawPen, rPSIMD);

        //Calculate [5]EPH  [8]EPV
        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.ePenH, 
            _mm256_mul_pd(xcos(angleRSIMD),rPSIMD)
        );

        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.ePenD, 
            _mm256_mul_pd(xcos(angleSIMD),rPSIMD)
        );
        
        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.ePenHN,
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleRSIMD, normalizationR))
            )
        );
        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.ePenDN,
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleSIMD, normalizationR))
            )
        );

        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.tToTarget,
            tSIMD
        );

        _mm256_storeu_pd(s.stdData.data()+i+s.sizeAligned*s.tToTargetA,
            _mm256_div_pd(
                tSIMD,
                _mm256_set1_pd(3.1)                 
            )
        );
    }
    

    public:

    double calcNormalizationR(const double angle, const double normalizationR){ //Input in radians
        if(fabs(angle) > normalizationR){
            return fabs(angle) - normalizationR;
        }else{
            return 0;
        }
    }

    #ifdef USE_SIMD

    __m256d calcNormalizationRSIMD(__m256d angle, const double normalizationR){
        return _mm256_max_pd(
            _mm256_sub_pd(
                _mm256_and_pd(
                    angle, _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF))
                ),
                _mm256_set1_pd(normalizationR)
            )
            , _mm256_set1_pd(0.0)
        );
    }

    #endif

    inline int signum(double x){
        return ((0.0) < x) - (x < (0.0));
    }

    shellCalc() = default;

    void editTestParameters(double max, double min, double precision, double x0, double y0, double dt){
        if(!max){
            this->max = max;
        }
        if(!min){
            this->min = max;
        }
        if(!precision){
            this->precision = max;
        }
        if(!x0){
            this->x0 = max;
        }
        if(!y0){
            this->y0 = max;
        }
        if(!dt){
            this->dt = max;
        }
    }

    void calculateStd(shell& s){
        unsigned int i;
        s.size = (unsigned int) (max - min) / precision;
        s.sizeAligned = (sizeof(__m256d)/sizeof(double) - (s.size % (sizeof(__m256d)/sizeof(double)))) + s.size;

        //oneVector.resize(size);
        //memset(oneVector.data(), 1, sizeof(double) * size);

        s.trajectories.resize(2 * s.size);
        s.stdData.resize(s.maxColumns * s.sizeAligned);

        //omp_set_num_threads(6);
        #pragma omp parallel for schedule(dynamic)
        for(i=0; i<s.size; i+=vSize){
            singleTraj(i, s, true);
            //printf("i %d \n", i);
        }
        s.completed = true;
    }

    //Post-Penetration 
    private:
    double threshold, fuseTime, dtf = 0.0001;
    double xf0 = 0, yf0 = 0;
    bool completed = false;

    void postPenTraj(const unsigned int i, shell& s, double v_x, double v_y, double v_z, double thickness){
        const procCharacteristic pC = s.retPC();
        double k = pC.k;
        double cw_2 = pC.cw_2;
        double pPPC = pC.pPPC;
        double normalizationR = pC.normalizationR;
        double T, p, rho, t, x, y, z; //v_x, v_y, v_z;
        //v_x = velocities[i];
        //v_z = velocities[i+postPenSizeAligned*2];
        //v_y = velocities[i+postPenSizeAligned];
        x = xf0;
        y = yf0;
        z = xf0;
        t = 0;
        if(v_x > 0){
            while(t < fuseTime){
                x = x + dtf*v_x;
                z = z + dtf*v_z;
                y = y + dtf*v_y;
                T = t0 - L*y;
                p = p0*pow((1-L*y/t0),(g*M/(R*L)));
                rho = p*M/(R*T);

                v_x = v_x - dtf*k*rho*(cw_1*v_x*v_x+cw_2*v_x);
                v_z = v_z - dtf*k*rho*(cw_1*v_z*v_z+cw_2*v_z);
                v_y = v_y - dtf*(g - k*rho*(cw_1*v_y*v_y+cw_2*fabs(v_y))*signum(v_y));
                t = t + dtf;
            }
            s.postPenData[i+s.postPenSize*2] = x;
            s.postPenData[i+s.postPenSize*3] = y;
            s.postPenData[i+s.postPenSize*4] = z;
            if(thickness > pC.threshold){
                s.postPenData[i+s.postPenSize*5] = x;
            }else{
                s.postPenData[i+s.postPenSize*5] = -1;
            }
        }else{
            s.postPenData[i+s.postPenSize*2] = 0;
            s.postPenData[i+s.postPenSize*3] = 0;
            s.postPenData[i+s.postPenSize*4] = 0;
            s.postPenData[i+s.postPenSize*5] = 0;
        }
    }

    public:
    bool includeNormalization = true;
    bool nChangeTrajectory = true;
    unsigned int postPenSize, postPenSizeAligned, angleSize;
    //std::vector<double> angles;

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE SIMD ALIGNMENT MAY NOT BE GUARANTEED
     * [0:1) Lateral Angle [1:2) Distance [2:3) X [3:4) Y [4:5) Z [5:6) XWF
     */
    //std::vector<double, AlignmentAllocator<double, 32>> postPenData;


    void calculatePostPen(const double thickness, shell& s){
        if(!s.completed){
            calculateStd(s);
        }

        unsigned int i, distIndex, anglesIndex;
        s.postPenSize = s.size * s.angles.size();
        s.postPenSizeAligned = s.sizeAligned * s.angles.size();
        s.postPenData.resize(6 * postPenSize);


        #ifdef USE_SIMD
        __m256d hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV, v0V, pPVV, ePenetrationV, eThickness, hFAngleV, vFAngleV, v_x, v_y, v_z;
        #endif
        /*for(unsigned int j=0; j<int64_vec::size; j++){
            temp.set_elt(j, j);
        }*/

        //#ifdef USE_SIMD
        #pragma omp parallel for private(distIndex, hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV, pPVV, ePenetrationV, hFAngleV, vFAngleV)
        for(unsigned int i=0; i < s.postPenSize % vSize; i+=vSize){
            distIndex = (i < s.size) ? i : i % s.size;
            anglesIndex = i / s.size;

            //hAngleV = postPenData.data() + i;
            unsigned int j, k = 0;
            if(distIndex < s.size - vSize + 1){
                hAngleV = _mm256_loadu_pd(s.angles.data()+anglesIndex);
                vAngleV = _mm256_loadu_pd(s.stdData.data()+s.sizeAligned*s.impactAHR);
                ePenetrationV = _mm256_loadu_pd(s.stdData.data()+s.sizeAligned*s.rawPen);
                v0V = _mm256_loadu_pd(s.stdData.data()+s.sizeAligned*s.impactV);
            }else{
                for(j = 0; (j + distIndex < s.size) && (j < vSize); j++){
                    hAngleV[j] = s.angles[anglesIndex];
                    vAngleV[j] = s.stdData[distIndex+j+s.sizeAligned*s.impactAHR];
                    ePenetrationV[j] = s.stdData[distIndex+j+s.sizeAligned*s.rawPen];
                    v0V[j] = s.stdData[distIndex+j+s.sizeAligned*s.impactV];
                }
                if(anglesIndex < s.angles.size()){
                    for(; (j < vSize); j++){
                        hAngleV[j] = s.angles[anglesIndex + 1];
                        vAngleV[k] = s.stdData[k+s.sizeAligned*s.impactAHR];
                        ePenetrationV[k] = s.stdData[k+s.sizeAligned*s.rawPen];
                        v0V[j] = s.stdData[k+s.sizeAligned*s.impactV];
                        k++;
                    }
                }
            }

            
            cAngleV = xacos(_mm256_mul_pd(xcos(hAngleV), xcos(vAngleV)));
            //cAngleV = acos(cos(hAngleV) * cos(vAngleV));
            nCAngleV = calcNormalizationRSIMD(cAngleV, s.retPC().normalizationR);
            eThickness = _mm256_div_pd(_mm256_set1_pd(thickness), xcos(nCAngleV));
                
            //ePenetrationV = _mm256_mul_pd(xcos(nCAngleV), ePenentrationV);
            pPVV = _mm256_max_pd(_mm256_mul_pd(v0V, 
                _mm256_sub_pd(_mm256_set1_pd(1), 
                    xexp(_mm256_sub_pd(_mm256_set1_pd(1),
                        _mm256_div_pd(ePenetrationV, eThickness)
                    )))
            ), _mm256_set1_pd(0));
            //pPVV = ifthen(ePenetrationV > float64_vec(thickness), 
            //(float64_vec(1)-exp(float64_vec(1)-ePenetrationV/float64_vec(thickness)))*v0V, float64_vec(0.0));

            //aAngleV = nCAngleV / cAngleV;
            //hFAngleV = hAngleV * aAngleV;
            //vFAngleV = vAngleV * aAngleV;

            aAngleV = _mm256_div_pd(nCAngleV, cAngleV);
            hFAngleV = _mm256_mul_pd(hAngleV, aAngleV);
            vFAngleV = _mm256_mul_pd(vAngleV, aAngleV);

            v_x = _mm256_mul_pd(pPVV,
                _mm256_mul_pd(xcos(vFAngleV), xcos(hFAngleV))
            );
            v_y = _mm256_mul_pd(pPVV, xsin(vFAngleV));
            v_z = _mm256_mul_pd(pPVV,
                _mm256_mul_pd(xcos(vFAngleV), xsin(hFAngleV))
            );

            for(j=0; (j<vSize) && (j+i < s.postPenSize); j++){
                postPenTraj(i, s, v_x[j], v_y[j], v_z[j], eThickness[j]);
            }

            //storea(pPVV * cos(vFAngleV) * cos(vFAngleV)  , velocities.data()+i);
            //storea(pPVV * sin(vFAngleV)                  , velocities.data()+i+  postPenSizeAligned);
            //storea(pPVV * cos(vFAngleV) * sin(hFAngleV)  , velocities.data()+i+2*postPenSizeAligned);
            //storea(float64_vec(thickness) / cos(nCAngleV), velocities.data()+i+3*postPenSizeAligned);
        }

    }


};

