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

#ifdef USE_SIMD
//#include "vecmathlib/vecmathlib.h"
#define ENABLE_AVX2
#include "sleef-2.80/simd/sleefsimddp.c"
//#include "alignmentAllocator.h"
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
    std::string name; 
}shipParams;

//using namespace vecmathlib;

class shell{
    private:
    double C = 0.5561613;
    double g = 9.81;
    double t0 = 288;
    double L = 0.0065;
    double p0 = 101325;
    double R = 8.31447;
    double M = 0.0289644;
    double cw_1 = 1;

    double v0; // = 780;
    double caliber; // = .460;
    double krupp; // = 2574;
    double mass; // = 1460;
    double normalization; // = 6;
    double cD; // = .292;
    std::string name;

    double max = 25;
    double min = 0;
    double precision = .1;
    double x0 = 0, y0 = 0;
    double dt = .01;

    double k, cw_2, pPPC, normalizationR;

    std::vector<double> oneVector;
    std::vector<double> temp;

    unsigned int alignmentRequired = vSize * 8;
    /*standardOut1 - indicies refer to multiples of size
    [0:1) v_x [1:2) v_y
    */

    void preProcess(){
        k = 0.5 * cD * pow((caliber/2),2) * M_PI / mass;
        cw_2 = 100+1000/3*caliber;
        pPPC = C * krupp/2400 * pow(mass,0.55) / pow((caliber*1000),0.65);
        normalizationR = normalization / 180 * M_PI;
    }

    //std::vector<double> stdOut0;

    void singleTraj(unsigned int i, bool addTraj){
        //std::cout<<"Running0 "<< i<<std::endl;
        __m256d angleSIMD, angleRSIMD, temp, v0SIMD = _mm256_set1_pd(v0);
        __m256d vx, vy, tSIMD;
        temp = {0.0, 1.0, 2.0, 3.0};

        angleSIMD = _mm256_fmadd_pd(
            _mm256_add_pd(_mm256_set1_pd(i), temp), 
            _mm256_set1_pd(precision),
            _mm256_set1_pd(min));
        angleRSIMD = _mm256_mul_pd(angleSIMD, _mm256_set1_pd((M_PI / 180)));
        _mm256_storeu_pd(stdData.data() + i + launchA * sizeAligned, angleSIMD);

        vx = _mm256_mul_pd(v0SIMD, xcos(angleRSIMD));
        vy = _mm256_mul_pd(v0SIMD, xsin(angleRSIMD));

        double T, p, rho, t, x, y, v_x, v_y;
        unsigned int counter;
        #define __TrajBuffer__ 128
        double xT[__TrajBuffer__], yT[__TrajBuffer__];

        for(unsigned int j = 0; (j+i<size) && (j < vSize); j++){

            if(addTraj){
                trajectories[2*(i+j)  ].reserve(__TrajBuffer__);
                trajectories[2*(i+j)+1].reserve(__TrajBuffer__);
            }
            //double T, p, rho, t, x, y, v_x, v_y;
            v_x = vx[j];
            v_y = vy[j];
            //std::cout<<vx[j]<<" "<<vy[j]<<std::endl;
            //std::cout<<k<<" "<<v_x<<" "<<v_y<<std::endl;
            x = x0;
            y = y0;
            t = 0;
            trajectories[2*(i+j)  ].push_back(x);
            trajectories[2*(i+j)+1].push_back(y);
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
                    trajectories[2*(i+j)  ].insert(trajectories[2*(i+j)  ].end(), xT, &xT[counter]);
                    trajectories[2*(i+j)+1].insert(trajectories[2*(i+j)+1].end(), yT, &yT[counter]);

                }
            }
            stdData[i+j+sizeAligned*distance] = x;
            //stdData[i+sizeAligned*tToTarget] = t;
            //stdOut0[i] = v_x;
            //stdOut0[i+sizeAligned] = v_y;
            vx[j] = v_x;
            vy[j] = v_y;
            tSIMD[j] = t;
            //std::cout<<k<<" "<<v_x<<std::endl;
        }

        /*std::cout<<vx[0]<<std::endl;
        std::cout<<vx[1]<<std::endl;
        std::cout<<vx[2]<<std::endl;
        std::cout<<vx[3]<<std::endl;*/
        //Calculate [2]IA , [7]IA_D
        __m256d iVSIMD, rPSIMD;  
        angleRSIMD = xatan(_mm256_div_pd(vy,vx));

        _mm256_storeu_pd(stdData.data()+i+sizeAligned*impactAHR, angleRSIMD);
        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*impactAHD, 
            _mm256_mul_pd(angleRSIMD, _mm256_set1_pd(180 / M_PI))
        );

        //iADRSIMD = float64_vec(M_PI / 2) + iARSIMD;
        angleSIMD = _mm256_add_pd(
            _mm256_set1_pd(M_PI / 2), angleRSIMD
        );

        _mm256_storeu_pd(stdData.data()+i+sizeAligned*impactADD,
            _mm256_mul_pd(angleSIMD, _mm256_set1_pd(180 / M_PI))
        );

        //Calculate [3]iV,  [4]rP

        iVSIMD = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(vx, vx), _mm256_mul_pd(vy, vy)));
        _mm256_storeu_pd(stdData.data()+i+sizeAligned*impactV, iVSIMD);

        rPSIMD = _mm256_mul_pd(xpow(iVSIMD, _mm256_set1_pd(1.1)), _mm256_set1_pd(pPPC));
        _mm256_storeu_pd(stdData.data()+i+sizeAligned*rawPen, rPSIMD);

        //Calculate [5]EPH  [8]EPV
        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*ePenH, 
            _mm256_mul_pd(xcos(angleRSIMD),rPSIMD)
        );

        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*ePenD, 
            _mm256_mul_pd(xcos(angleSIMD),rPSIMD)
        );
        
        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*ePenHN,
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleRSIMD))
            )
        );
        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*ePenDN,
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleSIMD))
            )
        );

        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*tToTarget,
            tSIMD
        );

        _mm256_storeu_pd(
            stdData.data()+i+sizeAligned*tToTargetA,
            _mm256_div_pd(
                tSIMD,
                _mm256_set1_pd(3.1)                 
            )
        );
    }
    

    public:

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

    //For convenience purposes
    /*
    const unordered_map<string, unsigned int> stdDataIndex = {
        {"distance"   ,  0}, {"launchA",  1}, {"impactAHR",  2},
        {"impactAHD"  ,  3}, {"impactV",  4}, {"rawPen"   ,  5},
        {"ePenH"      ,  6}, {"ePenHN" ,  7}, {"impactADD",  8},
        {"ePenD"      ,  9}, {"ePenDN" , 10}, {"tToTarget", 11},
        {"tToTarget-A", 12}
    }; */

    enum stdDataIndex{
        distance   , launchA, impactAHR, 
        impactAHD  , impactV, rawPen   , 
        ePenH      , ePenHN , impactADD, 
        ePenD      , ePenDN , tToTarget, 
        tToTargetA};

    std::vector<unsigned int> stdDataSizeIndex;

    //unordered_map<string, unsigned int> stdDataSizeIndex;

    unsigned int size;
    unsigned int sizeAligned; 

    double calcNormalizationR(const double angle){ //Input in radians
        if(fabs(angle) > normalizationR){
            return fabs(angle) - normalizationR;
        }else{
            return 0;
        }
    }

    #ifdef USE_SIMD
    /*
    float64_vec calcNormalizationRSIMD(const float64_vec angle){
        float64_vec result;
        std::vector<double, AlignmentAllocator<double,  vSize * 8>> aR;
        aR.resize(4);
        storea(angle, aR.data());
        double a;
        #pragma omp simd 
        for(int i=0; i<vSize; i++){
            a = aR[i];

            if(fabs(a) > normalizationR){
                result.set_elt(i, fabs(a) - normalizationR);
            }else{
                result.set_elt(i, 0);
            }
        }
        return result;
    }
    */
    /*float64_vec calcNormalizationRSIMD(float64_vec angle){
        return fmax(fabs(angle) - float64_vec(normalizationR), float64_vec(0.0));
    }*/

    __m256d calcNormalizationRSIMD(__m256d angle){
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

    /*
    template <typename T> inline constexpr
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }

    template <typename T> inline constexpr
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }

    template <typename T> inline constexpr
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }*/
    inline int signum(double x){
        return ((0.0) < x) - (x < (0.0));
    }

    shell();

    shell(shipParams sp){
        v0 = sp.v0;
        caliber = sp.caliber;
        krupp = sp.krupp;
        mass = sp.mass;
        normalization = sp.normalization;
        cD = sp.cD;
        name = sp.cD;

        if(sp.fuseTime){
            fuseTime = sp.fuseTime;
        }else{
            fuseTime = .033; 
        }

        if(sp.threshold){
            threshold = sp.threshold;
        }else{
            threshold = caliber / 6;
        }

        preProcess();
    }

    shell(double v0, double caliber, double krupp, double mass,
    double normalization, double cD, std::string name, double threshold, double fuseTime){
        this->fuseTime = fuseTime;

        this->v0 = v0;
        this->caliber = caliber;
        this->krupp = krupp;
        this->mass = mass;
        this->normalization = normalization;
        this->cD = cD;
        this->name = name;

        if(threshold){
            this->threshold = threshold;
        }else{
            this->threshold = caliber / 6;
        }

        preProcess();
        //assignShellValues(v0, caliber, krupp, mass, normalization, cD, name);
    }

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

    shipParams returnShipParams(){
        shipParams ret;
        ret.caliber = caliber;
        ret.cD = cD;
        ret.fuseTime = fuseTime;
        ret.krupp = krupp;
        ret.mass = mass;
        ret.name = name;
        ret.normalization = normalization;
        ret.threshold = threshold;
        ret.v0 = v0;
        return ret;
    }

    void calculateStd(){
        unsigned int i;
        size = (unsigned int) (max - min) / precision;
        sizeAligned = (sizeof(__m256d)/sizeof(double) - (size % (sizeof(__m256d)/sizeof(double)))) + size;

        oneVector.resize(size);
        memset(oneVector.data(), 1, sizeof(double) * size);

        trajectories.resize(2 * size);
        stdData.resize(13 * sizeAligned);

        //stdOut0.resize(sizeAligned * 2);

        stdDataSizeIndex.resize(13);
        for(i=0; i<13; i++){
            stdDataSizeIndex[i] = i * sizeAligned;
        }
        
        //omp_set_num_threads(6);
        #pragma omp parallel for //schedule(dynamic)
        for(i=0; i<size; i+=4){
            singleTraj(i, true);
            //printf("i %d \n", i);
        }


    }
    /*
    void printStdOut(){
        for(unsigned int i=0; i<size; i++){
            printf("%f %f %f %f %f\n", stdOut0[i], stdOut0[i+sizeAligned], 
            stdOut0[i+sizeAligned*2], stdOut0[i+sizeAligned*3], stdOut0[i+sizeAligned*4]);
        }
    }*/
    void printStdData(){
        for(unsigned int i=0; i<size; i++){
            for(unsigned int j=0; j<stdDataSizeIndex.size(); j++){
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
    /*
    double *exposeStdData(){
        return stdData.data();
    }*/

    std::vector<double> stdDataCopy(){
        std::vector<double> out;
        out.resize(size * stdDataSizeIndex.size());
        
        #pragma omp parallel for schedule(static)
        for(unsigned int i=0; i<stdDataSizeIndex.size(); i++){
            std::copy_n(stdData.begin() + i * sizeAligned, size, out.begin() + i * size);
        }

        //std::copy_n(stdData.begin(), stdData.size(), out.begin());
        //std::cout<<&out<<std::endl;
        return out;
    }


    //Post-Penetration 
    private:
    double threshold, fuseTime, dtf = 0.0001;
    double xf0 = 0, yf0 = 0;
    bool completed = false;
    //std::vector<double, AlignmentAllocator<double, 32>> velocities;
    
    /*
    void printVelocities(){
        printf("Velocities\n");
        for(unsigned int i=0; i<postPenSize; i++){
            for(int j=0; j<4; j++){
                printf("%f ", velocities[i + j * postPenSize]);
            }
            printf("\n");
        }
    }*/
#ifdef NOTTEST
    void postPenTraj(const unsigned int i, std::vector<double, AlignmentAllocator<double, 32>> velocities){
        double T, p, rho, t, x, y, z, v_x, v_y, v_z;
        v_x = velocities[i];
        v_z = velocities[i+postPenSizeAligned*2];
        v_y = velocities[i+postPenSizeAligned];
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
            postPenData[i+postPenSize*2] = x;
            postPenData[i+postPenSize*3] = y;
            postPenData[i+postPenSize*4] = z;
            if(velocities[i+postPenSizeAligned*3] > threshold){
                postPenData[i+postPenSize*5] = x;
            }else{
                postPenData[i+postPenSize*5] = -1;
            }
        }else{
            postPenData[i+postPenSize*2] = 0;
            postPenData[i+postPenSize*3] = 0;
            postPenData[i+postPenSize*4] = 0;
            postPenData[i+postPenSize*5] = 0;
        }
    }

    public:
    bool includeNormalization = true;
    bool nChangeTrajectory = true;
    unsigned int postPenSize, postPenSizeAligned, angleSize;
    std::vector<double> angles;

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE SIMD ALIGNMENT MAY NOT BE GUARANTEED
     * [0:1) Lateral Angle [1:2) Distance [2:3) X [3:4) Y [4:5) Z [5:6) XWF
     */
    std::vector<double, AlignmentAllocator<double, 32>> postPenData;


    void setAngles(std::vector<double> *anglesIn){
        //this->angles.reset(angles);
        //angleSize = this->angles->size();
        angleSize = anglesIn->size();
        angles.resize(angleSize);
        std::copy_n(anglesIn->begin(), angleSize, angles.begin());
    }

    void setAngles(std::vector<double> anglesIn){
        /*
        std::cout<<&anglesIn<<std::endl;
        this->angles.reset(&anglesIn);
        //this->angles.reset(new vector<>)
        //std::copy_n(angles, angles.size(), )
        std::cout<<this->angles<<std::endl;
        std::cout<<this->angles->size()<<std::endl;
        angleSize = this->angles->size();
        for(int i=0; i<angleSize; i++){
            std::cout<<this->angles->at(i)<<std::endl;
        }*/
        angleSize = anglesIn.size();
        angles.resize(angleSize);
        std::copy_n(anglesIn.begin(), angleSize, angles.begin());
    }

    void setAngles(double *anglesIn, unsigned int size){
        /*
        this->angles.reset(new std::vector<double>);
        this->angles->resize(size);
        std::copy_n(angles, size, this->angles->data());
        angleSize = this->angles->size();*/
        angleSize = size;
        angles.resize(size);
        std::copy_n(anglesIn, size, angles.begin());
    }

    void editPostPen(double dtf){
        if(dtf){
            this->dtf = dtf;
        }
    }

    void calculatePostPen(const double thickness){
        if(!completed){
            calculateStd();
        }

        unsigned int i;
        postPenSize = size * angleSize;
        postPenSizeAligned = sizeAligned * angleSize;
        postPenData.resize(6 * postPenSize);

        std::vector<double, AlignmentAllocator<double, 32>> velocities;
        velocities.resize(4 * postPenSizeAligned);

        //printPostPen();
        omp_set_num_threads(6);
        #pragma omp parallel for
        for(unsigned int i=0; i<angleSize; i++){
            //std::cout<<i<<std::endl;
            //std::cout<<angles[i]<<std::endl;
            fill(postPenData.begin() + i * size, postPenData.begin() + (i+1) * size, angles[i]);
            copy_n(stdData.begin() + sizeAligned * distance, size, postPenData.begin() + postPenSize + i * size);
        }

        double hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle;
        unsigned int distIndex, dAIT;
        //float64_vec distIndexV, anglesIndexV;
        #ifdef USE_SIMD
        float64_vec hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV, v0V, pPVV, ePenetrationV, hFAngleV, vFAngleV;
        #endif
        /*for(unsigned int j=0; j<int64_vec::size; j++){
            temp.set_elt(j, j);
        }*/


        if(includeNormalization){
            if(nChangeTrajectory){
                i=0;
                #ifdef USE_SIMD
                #pragma omp parallel for private(distIndex, hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV, pPVV, ePenetrationV, hFAngleV, vFAngleV)
                for(unsigned int i=0; i < (size * angleSize) - (size * angleSize) % vSize; i+=vSize){
                    distIndex = i % size;
                    //anglesIndex = i / size;

                    hAngleV = postPenData.data() + i;

                    if(size - distIndex >= vSize){
                        vAngleV =       stdData.data()+distIndex+sizeAligned*impactAHR;
                        ePenetrationV = stdData.data()+distIndex+sizeAligned*rawPen;
                        v0V =           stdData.data()+distIndex+sizeAligned*impactV;
                    }else{
                        for(unsigned int j=0; j<vSize; j++){
                            dAIT = (i + j) % size;
                            vAngleV.set_elt(      j, stdData[dAIT+sizeAligned*impactAHR]);
                            ePenetrationV.set_elt(j, stdData[dAIT+sizeAligned*rawPen]);
                            v0V.set_elt(          j, stdData[dAIT+sizeAligned*impactV]);
                        }
                    }

                    cAngleV = acos(cos(hAngleV) * cos(vAngleV));
                    nCAngleV = calcNormalizationRSIMD(cAngleV);
                        
                    ePenetrationV *= cos(nCAngleV);
                    pPVV = ifthen(ePenetrationV > float64_vec(thickness), 
                    (float64_vec(1)-exp(float64_vec(1)-ePenetrationV/float64_vec(thickness)))*v0V, float64_vec(0.0));

                    aAngleV = nCAngleV / cAngleV;
                    hFAngleV = hAngleV * aAngleV;
                    vFAngleV = vAngleV * aAngleV;

                    storea(pPVV * cos(vFAngleV) * cos(vFAngleV)  , velocities.data()+i);
                    storea(pPVV * sin(vFAngleV)                  , velocities.data()+i+  postPenSizeAligned);
                    storea(pPVV * cos(vFAngleV) * sin(hFAngleV)  , velocities.data()+i+2*postPenSizeAligned);
                    storea(float64_vec(thickness) / cos(nCAngleV), velocities.data()+i+3*postPenSizeAligned);
                }
                #else
                #pragma omp parallel for private(distIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
                #endif
                for(; i < size * angleSize; i++){
                    distIndex = i % size;
                    //anglesIndex = i / size;

                    hAngle = postPenData[i] /180*M_PI;
                    vAngle = stdData[distIndex+sizeAligned*impactAHR];
                    cAngle = acos(cos(hAngle) * cos(vAngle));
                    nCAngle = calcNormalizationR(cAngle);
 
                    ePenetration = stdData[distIndex+sizeAligned*rawPen]*cos(nCAngle);

                    if(ePenetration > thickness){
                        pPV = (1-exp(1-ePenetration/thickness)) * stdData[distIndex+sizeAligned*impactV];
                    }else{
                        pPV = 0;
                    }

                    aAngle = nCAngle / cAngle;
                    hFAngle = hAngle * aAngle;
                    vFAngle = vAngle * aAngle;

                    velocities[i                     ] = pPV * cos(vFAngle) * cos(vFAngle);
                    velocities[i+  postPenSizeAligned] = pPV * sin(vFAngle);
                    velocities[i+2*postPenSizeAligned] = pPV * cos(vFAngle) * sin(hFAngle);
                    velocities[i+3*postPenSizeAligned] = thickness / cos(nCAngle);
                }
            }else{
                i=0;
                #pragma omp parallel for private(distIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
                for(unsigned int i=0; i < (size * angleSize) - (size * angleSize) % vSize; i+=vSize){
                    distIndex = i % size;
                    //anglesIndex = i / size;

                    hAngleV = postPenData.data() + i;

                    if(size - distIndex >= vSize){
                        vAngleV = stdData.data()+distIndex+sizeAligned*impactAHR;
                        ePenetrationV = stdData.data()+distIndex+sizeAligned*rawPen;
                        v0V = stdData.data()+distIndex+sizeAligned*impactV;
                    }else{
                        for(unsigned int j=0; j<vSize; j++){
                            dAIT = (i + j) % size;
                            vAngleV.set_elt(j, stdData[dAIT+sizeAligned*impactAHR]);
                            ePenetrationV.set_elt(j, stdData[dAIT+sizeAligned*rawPen]);
                            v0V.set_elt(j, stdData[dAIT+sizeAligned*impactV]);
                        }
                    }

                    cAngleV = acos(cos(hAngleV) * cos(vAngleV));
                    nCAngleV = calcNormalizationRSIMD(cAngleV);
                        
                    ePenetrationV *= cos(nCAngleV);
                    pPVV = ifthen(ePenetrationV > float64_vec(thickness), (float64_vec(1)-exp(float64_vec(1)-ePenetrationV/float64_vec(thickness)))*v0V, float64_vec(0.0));

                    storea(pPVV * cos(vAngleV) * cos(vAngleV)    , velocities.data()+i);
                    storea(pPVV * sin(vAngleV)                   , velocities.data()+i+postPenSizeAligned);
                    storea(pPVV * cos(vAngleV) * sin(hAngleV)    , velocities.data()+i+2*postPenSizeAligned);
                    storea(float64_vec(thickness) / cos(nCAngleV), velocities.data()+i+3*postPenSizeAligned);
                }

                //#pragma omp parallel for private(distIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
                for(; i < size * angleSize; i++){
                    distIndex = i % size;
                    //anglesIndex = i / size;

                    hAngle = postPenData[i] /180*M_PI;
                    vAngle = stdData[distIndex+sizeAligned*impactAHR];
                    cAngle = acos(cos(hAngle) * cos(vAngle));
                    nCAngle = calcNormalizationR(cAngle);
                   
                    ePenetration = stdData[distIndex+sizeAligned*rawPen]*cos(nCAngle);

                    if(ePenetration > thickness){
                        pPV = (1-exp(1-ePenetration/thickness)) * stdData[distIndex+sizeAligned*impactV];
                    }else{
                        pPV = 0;
                    }

                    velocities[i                     ] = pPV * cos(vAngle) * cos(vAngle);
                    velocities[i+  postPenSizeAligned] = pPV * sin(vAngle);
                    velocities[i+2*postPenSizeAligned] = pPV * cos(vAngle) * sin(hAngle);
                    velocities[i+3*postPenSizeAligned] = thickness / cos(nCAngle);
                }
            }
        }else{
            i=0;
            #pragma omp parallel for private(distIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
            for(unsigned int i=0; i < (size * angleSize) - (size * angleSize) % vSize; i+=vSize){
                distIndex = i % size;
                //anglesIndex = i / size;

                hAngleV = postPenData.data() + i;

                if(size - distIndex >= vSize){
                    vAngleV = stdData.data()+distIndex+sizeAligned*impactAHR;
                    ePenetrationV = stdData.data()+distIndex+sizeAligned*rawPen;
                    v0V = stdData.data()+distIndex+sizeAligned*impactV;
                }else{
                    for(unsigned int j=0; j<vSize; j++){
                        dAIT = (i + j) % size;
                        vAngleV.set_elt(j, stdData[dAIT+sizeAligned*impactAHR]);
                        ePenetrationV.set_elt(j, stdData[dAIT+sizeAligned*rawPen]);
                        v0V.set_elt(j, stdData[dAIT+sizeAligned*impactV]);
                    }
                }

                cAngleV = acos(cos(hAngleV) * cos(vAngleV));
                //nCAngleV = calcNormalizationRSIMD(cAngle);
                    
                ePenetrationV *= cos(cAngleV);
                pPVV = ifthen(ePenetrationV > float64_vec(thickness), (float64_vec(1)-exp(float64_vec(1)-ePenetrationV/float64_vec(thickness)))*v0V, float64_vec(0.0));

                storea(pPVV * cos(vAngleV) * cos(vAngleV)  , velocities.data()+i);
                storea(pPVV * sin(vAngleV)                  , velocities.data()+i+postPenSizeAligned);
                storea(pPVV * cos(vAngleV) * sin(hAngleV)  , velocities.data()+i+2*postPenSizeAligned);
                storea(float64_vec(thickness) / cos(cAngleV), velocities.data()+i+3*postPenSizeAligned);
            }

            //#pragma omp parallel for private(distIndex, anglesIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
            for(; i < size * angleSize; i++){
                distIndex = i % size;
                //anglesIndex = i / size;

                hAngle = postPenData[i] /180*M_PI;
                vAngle = stdData[distIndex+sizeAligned*impactAHR];
                cAngle = acos(cos(hAngle) * cos(vAngle));
                //nCAngle = calcNormalizationR(cAngle);
                
                ePenetration = stdData[distIndex+sizeAligned*rawPen]*cos(cAngle);

                if(ePenetration > thickness){
                    pPV = (1-exp(1-ePenetration/thickness)) * stdData[distIndex+sizeAligned*impactV];
                }else{
                    pPV = 0;
                }

                velocities[i                     ] = pPV * cos(vAngle) * cos(vAngle);
                velocities[i+  postPenSizeAligned] = pPV * sin(vAngle);
                velocities[i+2*postPenSizeAligned] = pPV * cos(vAngle) * sin(hAngle);
                velocities[i+3*postPenSizeAligned] = thickness / cos(cAngle);
            }
        }

        #pragma omp parallel for
        for(unsigned int i=0; i<postPenSize; i++){
            postPenTraj(i, velocities);
        }
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

    std::vector<double> postPenDataCopy(){
        std::vector<double> out;
        out.resize(postPenData.size());
        std::copy_n(postPenData.begin(), postPenData.size(), out.begin());
        //std::cout<<&out<<std::endl;
        return out;
        //return postPenData;
    }
    /*
    double *exposePostPenData(){
        return postPenData.data();
    }

    double *exposePostPenDataCopy(){
        double *copy = new double[postPenData.size()];
        std::copy_n(postPenData.begin(), postPenData.size(), copy);
        return copy;
    }*/
#endif
};

