#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
//#include <omp.h>
//#include <immintrin.h>
//#include <x86intrin.h>
//#include "constant.h"

//#define USE_SIMD

#ifdef _SINGLE_PRECISION
typedef float fPType;
#else
typedef double fPType;
#endif

#include "concurrentqueue/concurrentqueue.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <cstring>
#include <algorithm>
#include <type_traits>
#include <thread>
#include <atomic>
#include <functional>
		
/*
double operator"" _kg (long double input){return input;}
double operator"" _lbs(long double input){return input * 0.453592;}

double operator"" _mps(long double input){return input;}
double operator"" _fps(long double input){return input * 0.3048;}
*/

/* Base shell characteristics 
 * May be used to implement a hash table in the future
 */
namespace shell{



namespace impact{
static constexpr unsigned int maxColumns = 13;
enum impactDataIndex{
    distance   , launchA, impactAHR, 
    impactAHD  , impactV, rawPen   , 
    ePenH      , ePenHN , impactADD, 
    ePenD      , ePenDN , tToTarget, 
    tToTargetA};
static_assert(tToTargetA == (maxColumns - 1), "Invalid standard columns");
}

namespace post{
static constexpr unsigned int maxColumns = 6;
enum postPenDataIndex{
    angle, distance, 
    x    , y       ,
    z    , xwf
};
static_assert(xwf == (maxColumns - 1), "Invaild postpen columns");
}

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

class shell{
    private:              //Description         units
    double v0;            //muzzle velocity     m/s
    double caliber;       //shell caliber       m
    double krupp;         //shell krupp         [ndim]
    double mass;          //shell mass          kg
    double cD;            //coefficient of drag [ndim]
    double normalization; //shell normalization degrees
    double threshold;     //fusing threshold    mm
    double fuseTime;      //fuse time           s
    std::string name;

    double k, cw_2, pPPC, normalizationR;

    //Condenses initial values into values used by calculations
    //[Reduces repeated computations]
    void preProcess(){
        k = 0.5 * cD * pow((caliber/2),2) * M_PI / mass; //condensed drag coefficient
        cw_2 = 100+1000/3*caliber; //quadratic drag coefficient
        pPPC = 0.5561613 * krupp/2400 * pow(mass,0.55) / pow((caliber*1000),0.65); //condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI; //normalization (radians)
    }

    public:
    unsigned int impactSize, postPenSize; //number of distances in: standard, postPen
    unsigned int impactSizeAligned, postPenSizeAligned; //Not 100% necessary - sizes adjusted to fulfill alignment
    bool completedImpact = false, completedPostPen = false;

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
    Refer to stdDataIndex enums defined above
    */
    std::vector<double> impactData;

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE SIMD ALIGNMENT MAY NOT BE GUARANTEED
     * [0:1) Lateral Angle [1:2) Distance 
     * [2:3) X             [3:4) Y 
     * [4:5) Z             [5:6) XWF
     * See enums defined above
     */
    std::vector<double> postPenData;
    
    shell() = default;

    shell(const double v0, const double caliber, const double krupp, const double mass,
    const double normalization, const double cD, const std::string& name, const double threshold, const double fuseTime = .033){
        setValues(v0, caliber, krupp, mass, normalization, cD, name, threshold, fuseTime);
    }

    void setValues(const double v0, const double caliber, const double krupp, const double mass,
    const double normalization, const double cD, const std::string& name, const double threshold, const double fuseTime = .033){
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
    }
    //Getter Functions
    double& getImpact(unsigned int i, unsigned int j){
        return impactData[i + j * impactSizeAligned];
    }

    double* getImpactPtr(unsigned int i, unsigned int j){
        return impactData.data() + i + j * impactSizeAligned;
    }

    double& getPostPen(unsigned int i, unsigned int j, unsigned int k){
        return postPenData[i + j * postPenSize + k * impactSize];
    }

    double* getPostPenPtr(unsigned int i, unsigned int j, unsigned int k){
        return postPenData.data() + i + j * postPenSize + k * impactSize;
    }

    const double& get_v0(){
        return v0;
    }
    const double& get_k(){
        return k;
    }
    const double& get_cw_2(){
        return cw_2;
    }
    const double& get_pPPC(){
        return pPPC;
    }
    const double& get_normalizationR(){
        return normalizationR;
    }
    const double& get_threshold(){
        return threshold;
    }
    const double& get_fuseTime(){
        return fuseTime;
    }
    
    //Could be split up into two classes
    shellParams returnShipParams(){
        shellParams ret;
        ret.caliber = caliber;
        ret.cD = cD;
        ret.fuseTime = fuseTime;
        ret.krupp = krupp;
        ret.mass = mass;
        //ret.name = name;
        ret.normalization = normalization;
        ret.threshold = threshold;
        ret.v0 = v0;
        return ret;
    }

    void printPostPenData(){
        for(unsigned int i=0; i<postPenSize; i++){
            for(unsigned int j=0; j<post::maxColumns; j++){
                std::cout<< std::fixed<< std::setprecision(4) << postPenData[i + j * postPenSize] << " "; 
            }
            std::cout<<"\n";
        }
        std::cout<<"Completed Post-Penetration\n";
    }

    void printImpactData(){
        for(unsigned int i=0; i<impactSize; i++){
            for(unsigned int j=0; j<impact::maxColumns; j++){
                std::cout<<std::fixed<<std::setprecision(4)<<getImpact(i, j)<< " ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"Completed Standard Data"<<std::endl;
    }
    void printTrajectory(unsigned int target){
        if(target >= impactSize){
            std::cout<<"Target Not Within Range of: "<<impactSize<< std::endl;
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
    //Threading
    std::atomic<int> counter, threadCount;
    int assigned, length;
    moodycamel::ConcurrentQueue<int> workQueue;
    static constexpr int workQueueBufferSize = 16;

    //Physical Constants     Description                    Units
    double g = 9.81;         //Gravitational Constant       m/(s^2)
    double t0 = 288;         //Temperature at Sea Level     K
    double L = 0.0065;       //Atmospheric Lapse Rate       C/m
    double p0 = 101325;      //Pressure at Sea Level        Pa
    double R = 8.31447;      //Ideal Gas Constant           J/(mol K)
    double M = 0.0289644;    //Molarity of Air at Sea Level kg/mol
    double cw_1 = 1;
    //Calculation Parameters 
    double max = 25;         //Max Angle                    degrees
    double min = 0;          //Min Angle                    degrees
    double precision = .1;   //Angle Step                   degrees
    double x0 = 0, y0 = 0;   //Starting x0, y0              m
    double dt = .01;         //Time step                    s

    //For vectorization - though not 100% necessary anymore since intrinsics were removed
    static_assert(sizeof(double) == 8, "Size of double is not 8 - required for AVX2"); //Use float64 in the future
    static_assert(std::numeric_limits<double>::is_iec559, "Type is not IEE754 compliant");
    static constexpr unsigned int vSize = (256 / 8) / sizeof(double);
    
    template<bool AddTraj>
    void singleTraj(const unsigned int i, const unsigned int j, shell&s, double* vx, double* vy, double* tVec){
        static constexpr unsigned int __TrajBuffer__ = 128;
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();

        double T, p, rho, t; //x, y, v_x, v_y;
        //__m128d pos, velocity, velocitySquared, dragIntermediary;
        double pos[2], velocity[2];
        int counter;
        if constexpr(AddTraj){
            s.trajectories[2*(i+j)  ].reserve(__TrajBuffer__);
            s.trajectories[2*(i+j)+1].reserve(__TrajBuffer__);
        }
        double xT[__TrajBuffer__], yT[__TrajBuffer__];

        //setting initial values
        velocity[0] = vx[j];                         //x component of velocity v_x
        velocity[1] = vy[j];                         //y component of velocity v_y
        pos[0] = x0;                                 //x start x0
        pos[1] = y0;                                 //y start y0
        s.trajectories[2*(i+j)  ].push_back(x0);     //add x start (x0) to trajectories
        s.trajectories[2*(i+j)+1].push_back(y0);     //add y start (y0) to trajectories
        t = 0;                                       //t start

        while(pos[1] >= 0){
            for(counter = 0; counter < __TrajBuffer__ && pos[1] >= 0; counter++){
                for(int l=0; l<2; l++){
                    pos[l] += velocity[l] * dt;
                }

                //Calculating air density
                T = t0 - L*pos[1];                       //Calculating air temperature at altitude
                p = p0*pow((1-L*pos[1]/t0),(g*M/(R*L))); //Calculating air pressure at altitude
                rho = p*M/(R*T);                         //Use ideal gas law to calculate air density

                //Calculate drag deceleration
                double dragIntermediary[2], velocitySquared[2];
                for(int l=0; l<2; l++){
                    velocitySquared[l] = velocity[l] * velocity[l]; //v^2 = v * v
                }

                dragIntermediary[0] = k*rho*(cw_1*velocitySquared[0] + cw_2*velocity[0]);                             //for horizontal (x) component
                dragIntermediary[1] = g - k*rho*(cw_1*velocitySquared[1]+cw_2*fabs(velocity[1]))*signum(velocity[1]); //for vertical   (y) component

                //Adjust for next cycle
                for(int l=0; l<2; l++){ //v -= drag * dt
                    velocity[l] -= dragIntermediary[l] * dt;
                }
                t += dt; //adjust time
                if constexpr(AddTraj){
                    xT[counter] = pos[0];                                                      
                    yT[counter] = pos[1];
                }
            }
            if constexpr(AddTraj){
                s.trajectories[2*(i+j)  ].insert(s.trajectories[2*(i+j)  ].end(), xT, &xT[counter]);
                s.trajectories[2*(i+j)+1].insert(s.trajectories[2*(i+j)+1].end(), yT, &yT[counter]);
            }

        }
        s.getImpact(i + j, impact::distance) = pos[0];
        vx[j] = velocity[0];
        vy[j] = velocity[1];
        tVec[j] = t;
    }

    template<bool AddTraj>
    void multiTraj(const unsigned int i, shell& s){
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_normalizationR();

        double vx[vSize], vy[vSize], tVec[vSize];
        for(int j=0; j<vSize; j++){
            s.getImpact(i + j, impact::launchA) = min + precision * (i+ j);

            double radianLaunch = s.getImpact(i + j, impact::launchA) * M_PI / 180;
            vx[j] = s.get_v0() * cos(radianLaunch);
            vy[j] = s.get_v0() * sin(radianLaunch);

        }

        for(int j = 0; (j+i<s.impactSize) && (j < vSize); j++){
            singleTraj<AddTraj>(i, j, s, vx, vy, tVec);
        }   

        for(int j=0; j<vSize; j++){
            double IA_R = atan(vy[j] / vx[j]);
            s.getImpact(i+j, impact::impactAHR) = IA_R;
            double IAD_R = M_PI / 2 + IA_R;
            double IA_D = IA_R * 180 / M_PI;
            s.getImpact(i+j, impact::impactAHD) = IA_D;
            s.getImpact(i+j, impact::impactADD) = 90 + IA_D;

            double IV = sqrt(vx[j] * vx[j] + vy[j] * vy[j]);
            s.getImpact(i+j, impact::impactV) = IV;
            double rawPen = pPPC * pow(IV, 1.1);
            s.getImpact(i+j, impact::rawPen) = rawPen;

            s.getImpact(i+j, impact::ePenH) = rawPen * cos(IA_R);
            s.getImpact(i+j, impact::ePenD) = rawPen * cos(IAD_R);

            s.getImpact(i+j, impact::ePenHN) = rawPen * cos(calcNormalizationR(IA_R, normalizationR));
            s.getImpact(i+j, impact::ePenDN) = rawPen * cos(calcNormalizationR(IAD_R, normalizationR));

            s.getImpact(i+j, impact::tToTarget) = tVec[j];
            s.getImpact(i+j, impact::tToTargetA) = tVec[j] / 3.1;
        }
    }
    
    template<bool AddTraj>
    void impactWorker(int threadId, shell *s){ //threadID is largely there for debugging
        while(counter < length){
            int index;
            if(workQueue.try_dequeue(index)){
                multiTraj<AddTraj>(index, *s);
                counter.fetch_add(1, std::memory_order_relaxed);
            }else{
                std::this_thread::yield();
            }
        }
        threadCount.fetch_add(1, std::memory_order_relaxed);
    }

    public:
    double calcNormalizationR(const double angle, const double normalizationR){ //Input in radians
        return (fabs(angle) > normalizationR) * (fabs(angle) - normalizationR);
    }

    inline int signum(double x){
        return ((0.0) < x) - (x < (0.0));
    }

    shellCalc() = default;
    //Replace with setter in the future
    void editTestParameters(double max, double min, double precision, double x0, double y0, double dt){
        if(!max){
            this->max = max;
        }
        if(!min){
            this->min = min;
        }
        if(!precision){
            this->precision = precision;
        }
        if(!x0){
            this->x0 = x0;
        }
        if(!y0){
            this->y0 = y0;
        }
        if(!dt){
            this->dt = dt;
        }
    }

    void calculateImpact(shell& s, bool addTraj, unsigned int nThreads=std::thread::hardware_concurrency()){
        if(addTraj){
            calculateImpact<true>(s, nThreads);
        }else{
            calculateImpact<false>(s, nThreads);
        }
    }

    template<bool AddTraj>
    void calculateImpact(shell& s, unsigned int nThreads=std::thread::hardware_concurrency()){
        s.impactSize = (unsigned int) (max - min) / precision;
        s.impactSizeAligned = vSize - (s.impactSize % vSize) + s.impactSize;
        s.trajectories.resize(2 * s.impactSize);
        s.impactData.resize(impact::maxColumns * s.impactSizeAligned);

        //replace with unified thread manager - so we don't make and delete threads everytime something is run
        if(nThreads  > std::thread::hardware_concurrency()){
            nThreads = std::thread::hardware_concurrency();
        }
        length = ceil((double) s.impactSize / vSize);
        if(length > nThreads){
            assigned = nThreads;
        }else{
            assigned = length;
        }

        //assigned = 1;
        counter = 0;
        threadCount = 0;
        std::vector<std::thread> threads;
        shell* sPtr = &s;
        threads.reserve(assigned - 1);

        for(int i=0; i<assigned - 1; i++){
            threads.emplace_back([=]{impactWorker<AddTraj>(i, sPtr);});
        }

        int buffer[workQueueBufferSize];
        int bCounter = 0;
        for(int i=0; i<s.impactSize; i+=vSize){
            buffer[bCounter] = i;
            bCounter++;
            if(bCounter == workQueueBufferSize){
                workQueue.enqueue_bulk(buffer, bCounter);
                bCounter = 0;
            }
            //multiTraj(i, s, addTraj);
        }
        workQueue.enqueue_bulk(buffer, bCounter);
        impactWorker<AddTraj>(assigned - 1, sPtr);
        

        while(threadCount < assigned){
            std::this_thread::yield();
        }
        //std::cout<<s.impactSize<<" "<<length<<" "<<counter<<"\n";
        for(int i=0; i<assigned-1; i++){
            threads[i].join();
        }
        s.completedImpact = true;
    }

    //Post-Penetration Section
    private:
    double dtf = 0.0001;
    double xf0 = 0, yf0 = 0;

    template<bool fast>
    void postPenTraj(const unsigned int i, shell& s, double v_x, double v_y, double v_z, double thickness){
        if constexpr(fast){
            double x = v_x * s.get_fuseTime();
            s.getPostPen(i, post::x  , 0) = x;
            s.getPostPen(i, post::y  , 0) = v_y * s.get_fuseTime();
            s.getPostPen(i, post::z  , 0) = v_z * s.get_fuseTime();
            s.getPostPen(i, post::xwf, 0) = (thickness >= s.get_threshold()) * x + !(thickness >= s.get_threshold()) * -1;
        }else{
            const double k = s.get_k();
            const double cw_2 = s.get_cw_2();
            double T, p, rho, t; 

            /* [indices]           0     1     2      
            * pos                 x     y     z      
            * velocities          v_x   v_y   v_z    
            * velocitiesSquared   v_x^2 v_y^2 v_z^2  
            * dragIntermediary    ad_x  ad_y  ad_z   
            * xz_dragIntermediary ad_x  ad_z       
            */
            double pos[3], velocities[3], velocitiesSquared[3], dragIntermediary[3];
            double xz_dragIntermediary[2]; 
            pos[0] = xf0, pos[1] = yf0, pos[2] = xf0;
            velocities[0] = v_x, velocities[1] = v_y, velocities[2] = v_z;
            t = 0;
            if(v_x > 0){
                while(t < s.get_fuseTime()){
                    for(int l=0; l<3; l++){
                        pos[l] += velocities[l] * dtf;
                    }
                    
                    //Calculate air density - likely unnecessary for this section as distances are so short
                    T = t0 - L*pos[1];
                    p = p0*pow((1-L*pos[1]/t0),(g*M/(R*L)));
                    rho = p*M/(R*T);

                    //Calculated drag deceleration
                    
                    for(int l=0; l<3; l++){
                        velocitiesSquared[l] = velocities[l] * velocities[l];
                    }
                    //velocitiesSquared = _mm256_mul_pd(velocities, velocities); //velocitiesSquared = velocities * velocities

                    xz_dragIntermediary[0] = (k * rho) * (cw_1 * velocitiesSquared[0] + cw_2 * velocities[0]);
                    xz_dragIntermediary[1] = (k * rho) * (cw_1 * velocitiesSquared[2] + cw_2 * velocities[2]);
                    //xz_dragIntermediary = (k * rho) * (cw_1 * velocitiesSquared[2, 0] + cw_2 * velocities[2, 0])
                    
                    dragIntermediary[0] = xz_dragIntermediary[0]; //x
                    dragIntermediary[1] = (g - k*rho*(cw_1*velocitiesSquared[1]+cw_2*fabs(velocities[1]))*signum(velocities[1])); 
                    dragIntermediary[2] = xz_dragIntermediary[1]; //z

                    //velocities -= dtf * dragIntermediary
                    for(int l=0; l<3; l++){
                        velocities[l] -= dtf * dragIntermediary[l];
                    }
                    t += dtf;                                                                             
                }
                s.getPostPen(i, post::x  , 0) = pos[0];
                s.getPostPen(i, post::y  , 0) = pos[1];
                s.getPostPen(i, post::z  , 0) = pos[2];
                s.getPostPen(i, post::xwf, 0) = (thickness >= s.get_threshold()) * pos[0] + !(thickness >= s.get_threshold()) * -1;
            }else{
                s.getPostPen(i, post::x  , 0) = 0;
                s.getPostPen(i, post::y  , 0) = 0;
                s.getPostPen(i, post::z  , 0) = 0;
                s.getPostPen(i, post::xwf, 0) = 0;
            }
        }
    }

    template<bool changeDirection, bool fast>
    void multiPostPen(int i, const double thickness, const double inclination_R, shell& s){
        //std::cout<<index<<"\n";
        //unsigned int i = index * vSize;
        double hAngleV[vSize], vAngleV[vSize];
        double v0V[vSize], penetrationV[vSize], eThicknessV[vSize];
        double v_x[vSize], v_y[vSize], v_z[vSize];
        unsigned int distIndex = (i < s.impactSize) ? i : i % s.impactSize;
        unsigned int anglesIndex = i / s.impactSize;

        unsigned int j, k = 0;

        if(i + vSize <= s.postPenSize){
            std::copy_n(s.getPostPenPtr(i, post::angle, 0), vSize, hAngleV);
        }else{
            for(j = 0; (i + j)< s.postPenSize; j++){
                hAngleV[j] = s.postPenData[i+j];
            }
        }

        if(distIndex < s.impactSize - vSize + 1){
            std::copy_n(s.getImpactPtr(distIndex, impact::impactAHR), vSize, vAngleV);
            std::copy_n(s.getImpactPtr(distIndex, impact::rawPen)   , vSize, penetrationV);
            std::copy_n(s.getImpactPtr(distIndex, impact::impactV)  , vSize, v0V);
        }else{
            for(j = 0; (j + distIndex < s.impactSize) && (j < vSize); j++){
                vAngleV[j] = s.getImpact(distIndex + j, impact::impactAHR);
                penetrationV[j] = s.getImpact(distIndex + j, impact::rawPen);
                v0V[j] = s.getImpact(distIndex + j, impact::impactV);
            }
            if(anglesIndex < s.postPenSize / s.impactSize){
                for(; (j < vSize); j++){
                    vAngleV[j] = s.getImpact(k, impact::impactAHR);
                    penetrationV[j] = s.getImpact(k, impact::rawPen);
                    v0V[j] = s.getImpact(k, impact::impactV);
                    k++;
                }
            }
        }

        for(int l=0; l<vSize; l++){
            double HA_R = hAngleV[l] * M_PI / 180;  //lateral  angle radians
            double VA_R = vAngleV[l] + inclination_R; //vertical angle radians
            double cAngle = acos(cos(HA_R) * cos(VA_R));
            double nCAngle = calcNormalizationR(cAngle, s.get_normalizationR());

            double eThickness = thickness / cos(nCAngle);
            double pPV = v0V[l] * (1 - exp(1 - penetrationV[l] / eThickness));

            if constexpr(changeDirection){
                double hFAngle = atan(tan(nCAngle) * tan(HA_R) / tan(cAngle));
                double vFAngle = atan(tan(nCAngle) * cos(hFAngle) * tan(VA_R) / cos(HA_R) / tan(cAngle));

                double v_x0 = pPV * cos(vFAngle) * cos(hFAngle);
                double v_y0 = pPV * cos(vFAngle) * sin(hFAngle);

                v_x[l] = v_x0 * cos(inclination_R) + v_y0 * sin(inclination_R);
                v_z[l] = v_y0 * cos(inclination_R) + v_x0 * sin(inclination_R);
                v_y[l] = pPV * sin(vFAngle);
            }else{
                v_x[l] = pPV * cos(VA_R) * cos(HA_R);
                v_z[l] = pPV * cos(VA_R) * sin(HA_R);
                v_y[l] = pPV * sin(VA_R);
            }

            eThicknessV[l] = eThickness;
        }

        for(unsigned int j=0; (j<vSize) && (j+i < s.postPenSize); j++){
            postPenTraj<fast>(i+j, s, v_x[j], v_y[j], v_z[j], eThicknessV[j]);
        }
        //std::cout<<index<<" Completed\n";
    }

    //Probably unnecessary... 
    void fillCopy(int id, shell* s, std::vector<double>* angles){
        //std::cout<<id<<"\n";
        for(int i=angles->size() * id / assigned; i<angles->size() * (id + 1) / assigned; i++){
            std::fill_n(s->postPenData.begin() + i * s->impactSize, s->impactSize, (double) angles->at(i));
            std::copy_n(s->getImpactPtr(0, impact::distance), s->impactSize, s->postPenData.begin() + s->postPenSize + i * s->impactSize);
        }
        counter.fetch_add(1, std::memory_order_relaxed);
    }

    void parallelFillCopy(shell* s, std::vector<double>* angles, unsigned int nThreads){
        std::vector<std::thread> threads;
        counter = 0;
        length = angles->size();
        if(length < nThreads){
            assigned = length;
        }else{
            assigned = nThreads;
        }
        for(int i=0; i<assigned - 1; i++){
            threads.push_back(std::thread([=]{fillCopy(i, s, angles);} ) );
        }
        fillCopy(assigned - 1, s, angles);
        while(counter < assigned){
            std::this_thread::yield();
        }
        for(int i=0; i<assigned - 1; i++){
            threads[i].join();
        }
    }

    public:
    bool includeNormalization = true;
    bool nChangeTrajectory = true;

    void calculatePostPen(const double thickness, const double inclination, shell& s, std::vector<double>& angles, const bool changeDirection=true,
    const bool fast=false, const unsigned int nThreads=std::thread::hardware_concurrency()){
        if(changeDirection){
            calculatePostPen<true>( thickness, inclination, s, angles, fast, nThreads);
        }else{
            calculatePostPen<false>(thickness, inclination, s, angles, fast, nThreads);
        }
    }
    
    template<bool changeDirection>
    void calculatePostPen(const double thickness, const double inclination, shell& s, std::vector<double>& angles, 
    const bool fast=false, const unsigned int nThreads=std::thread::hardware_concurrency()){
        if(fast){
            calculatePostPen<changeDirection, true>( thickness, inclination, s, angles, nThreads);
        }else{
            calculatePostPen<changeDirection, false>(thickness, inclination, s, angles, nThreads);
        }
    }

    template<bool changeDirection, bool fast>
    void calculatePostPen(const double thickness, const double inclination, shell& s, std::vector<double>& angles, const unsigned int nThreads){

        if(!s.completedImpact){
            std::cout<<"Standard Not Calculated - Running automatically\n";
            calculateImpact(s, false);
        }

        s.postPenSize = s.impactSize * angles.size();
        s.postPenData.resize(6 * s.postPenSize);

        parallelFillCopy(&s, &angles, nThreads);        

        std::vector<std::thread> threads;
        counter = 0, threadCount = 0;
        length = ceil((double) s.postPenSize/vSize);
        if(length < nThreads){
            assigned = length;
        }else{
            assigned = nThreads;
        }
        int buffer[workQueueBufferSize];
        int bCounter = 0;
        for(int i=0; i<s.postPenSize; i+=vSize){
            buffer[bCounter] = i;
            bCounter++;
            if(bCounter == workQueueBufferSize){
                workQueue.enqueue_bulk(buffer, bCounter);
                bCounter = 0;
            }
        }
        workQueue.enqueue_bulk(buffer, bCounter);
        double inclination_R = M_PI / 180 * inclination;
        shell* ptr = &s;
        for(int i=0; i<assigned - 1; i++){
            threads.emplace_back([=]{postPenWorker<changeDirection, fast>(i, thickness, inclination_R, ptr);});
        }
        postPenWorker<changeDirection, fast>(assigned - 1, thickness, inclination_R, &s);

        while(threadCount < assigned){
            std::this_thread::yield();
        }

        for(int i=0; i<assigned - 1; i++){
            threads[i].join();
        }
        s.completedPostPen = true;
        
    }

    private: 
    template<bool changeDirection, bool fast>
    void postPenWorker(int threadID, const double thickness, const double inclination, shell* s){
        while(counter < length){
            int index;
            if(workQueue.try_dequeue(index)){
                multiPostPen<changeDirection, fast>(index, thickness, inclination, *s);
                counter.fetch_add(1, std::memory_order_relaxed);
            }
            else{
                std::this_thread::yield();
            }
        }
        threadCount.fetch_add(1, std::memory_order_relaxed);
    }
};

}