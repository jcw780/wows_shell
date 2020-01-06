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
//const unsigned int vSize = 32 / sizeof(fPType);

/*
#ifdef _SINGLE_PRECISION
typedef __m256 fVType
#else
typedef __m256d fVType
#endif
*/

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

    double& getPostPen(unsigned int i, unsigned int j){
        return postPenData[i + j * postPenSize];
    }

    double* getPostPenPtr(unsigned int i, unsigned int j){
        return postPenData.data() + i + j * postPenSize;
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

    static_assert(sizeof(double) == 8, "Size of double is not 8 - required for AVX2"); //Use float64 in the future
    static_assert(std::numeric_limits<double>::is_iec559, "Type is not IEE754 compliant");
    static constexpr unsigned int vSize = (256 / 8) / sizeof(double);
    
    //template<bool AT>
    //void singleTraj(const unsigned int i, const unsigned int j, shell&s, __m256d& vx, __m256d& vy, __m256d& tSIMD);

    //template<>
    void singleTraj(const unsigned int i, const unsigned int j, shell&s, __m256d& vx, __m256d& vy, __m256d& tSIMD){
        static constexpr unsigned int __TrajBuffer__ = 128;
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();

        double T, p, rho, t; //x, y, v_x, v_y;
        __m128d pos, velocity, velocitySquared, dragIntermediary;
        unsigned int counter;
        s.trajectories[2*(i+j)  ].reserve(__TrajBuffer__);
        s.trajectories[2*(i+j)+1].reserve(__TrajBuffer__);
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
                pos = _mm_fmadd_pd(_mm_set1_pd(dt), velocity, pos); //positions += velocity * dt
                //Calculating air density
                T = t0 - L*pos[1];                       //Calculating air temperature at altitude
                p = p0*pow((1-L*pos[1]/t0),(g*M/(R*L))); //Calculating air pressure at altitude
                rho = p*M/(R*T);                         //Use ideal gas law to calculate air density

                //Calculate drag deceleration
                velocitySquared = _mm_mul_pd(velocity, velocity);                                                     //v^2 = v * v
                dragIntermediary[0] = k*rho*(cw_1*velocitySquared[0] + cw_2*velocity[0]);                             //for horizontal (x) component
                dragIntermediary[1] = g - k*rho*(cw_1*velocitySquared[1]+cw_2*fabs(velocity[1]))*signum(velocity[1]); //for vertical   (y) component

                //Adjust for next cycle
                velocity = _mm_fmadd_pd(_mm_set1_pd(-1 * dt), dragIntermediary, velocity); //v -= drag * dt
                t += dt;                                                                   //adjust time
                xT[counter] = pos[0];                                                      
                yT[counter] = pos[1];
            }
            s.trajectories[2*(i+j)  ].insert(s.trajectories[2*(i+j)  ].end(), xT, &xT[counter]);
            s.trajectories[2*(i+j)+1].insert(s.trajectories[2*(i+j)+1].end(), yT, &yT[counter]);

        }
        s.getImpact(i + j, impact::distance) = pos[0];
        vx[j] = velocity[0];
        vy[j] = velocity[1];
        tSIMD[j] = t;
    }

    void multiTraj(const unsigned int i, shell& s, const bool addTraj){
        //int i = index * vSize;
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_pPPC();
        __m256d angleSIMD, angleRSIMD, temp, v0SIMD = _mm256_set1_pd(s.get_v0());
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
        _mm256_storeu_pd(s.getImpactPtr(i, impact::launchA), angleSIMD);

        vx = _mm256_mul_pd(v0SIMD, xcos(angleRSIMD));
        vy = _mm256_mul_pd(v0SIMD, xsin(angleRSIMD));

        for(unsigned int j = 0; (j+i<s.impactSize) && (j < vSize); j++){
            singleTraj(i, j, s, vx, vy, tSIMD);
        }

        //Calculate [2]Impact Angles (impactAHR) [3]Impact Angles (impactAHD) [7] Impact Angle Deck (impactADD)
        __m256d iVSIMD, rPSIMD;  
        angleRSIMD = xatan(_mm256_div_pd(vy,vx));                 // = atan(vy / vx)

        _mm256_storeu_pd(s.getImpactPtr(i, impact::impactAHR) , angleRSIMD);
        _mm256_storeu_pd(s.getImpactPtr(i, impact::impactAHD),   
            _mm256_mul_pd(angleRSIMD, _mm256_set1_pd(180 / M_PI)) // = angleRSIMD * 180 / M_PI (convert to degrees) 
        );

        angleSIMD = _mm256_add_pd(_mm256_set1_pd(M_PI / 2), angleRSIMD);
        _mm256_storeu_pd(s.getImpactPtr(i, impact::impactADD),
            _mm256_mul_pd(angleSIMD, _mm256_set1_pd(180 / M_PI))
        );

        //Calculate [4]Impact Velocity (impactV),  [5]Raw Penetration (rawPen)
        iVSIMD = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(vx, vx), _mm256_mul_pd(vy, vy)));
        _mm256_storeu_pd(s.getImpactPtr(i, impact::impactV), iVSIMD);
        rPSIMD = _mm256_mul_pd(xpow(iVSIMD, _mm256_set1_pd(1.1)), _mm256_set1_pd(pPPC));
        _mm256_storeu_pd(s.getImpactPtr(i, impact::rawPen), rPSIMD);

        //Calculate [6]EPH  [9]EPV
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenH), _mm256_mul_pd(xcos(angleRSIMD),rPSIMD));
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenD), _mm256_mul_pd(xcos(angleSIMD),rPSIMD));
        
        //Calculate [7]EPHN [10]EPDN
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenHN),
            _mm256_mul_pd(rPSIMD, xcos(calcNormalizationRSIMD(angleRSIMD, normalizationR)))
        );
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenDN),
            _mm256_mul_pd(rPSIMD, xcos(calcNormalizationRSIMD(angleSIMD, normalizationR)))
        );

        //Calculate [11]TTT [12]TTTA
        _mm256_storeu_pd(s.getImpactPtr(i, impact::tToTarget), tSIMD);
        _mm256_storeu_pd(s.getImpactPtr(i, impact::tToTargetA), _mm256_div_pd(tSIMD, 
            _mm256_set1_pd(3.1))
        );
    }
    
    void impactWorker(int id, shell *s, bool addTraj){
        while(counter < length){
            int index;
            if(workQueue.try_dequeue(index)){
                multiTraj(index, *s, addTraj);
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

    #ifdef USE_SIMD

    inline __m256d calcNormalizationRSIMD(const __m256d angle, const double normalizationR){
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

    void calculateImpact(shell& s, bool addTraj){
        //unsigned int i;
        s.impactSize = (unsigned int) (max - min) / precision;
        s.impactSizeAligned = vSize - (s.impactSize % vSize) + s.impactSize;

        s.trajectories.resize(2 * s.impactSize);
        s.impactData.resize(impact::maxColumns * s.impactSizeAligned);

        //omp_set_num_threads(6);
        //#pragma omp parallel for schedule(dynamic, 2)
        length = s.impactSize / vSize;

        if(length > std::thread::hardware_concurrency()){
            assigned = std::thread::hardware_concurrency();
        }else{
            assigned = length;
        }
        counter = 0;
        threadCount = 0;
        std::vector<std::thread> threads;
        shell* sPtr = &s;
        for(int i=0; i<assigned - 1; i++){
            threads.emplace_back([=]{impactWorker(i, sPtr, addTraj);});
        }

        int buffer[32];
        int bCounter = 0;
        for(int i=0; i<s.impactSize; i+=vSize){
            buffer[bCounter] = i;
            bCounter++;
            if(bCounter == 32){
                workQueue.enqueue_bulk(buffer, bCounter);
                bCounter = 0;
            }
            //multiTraj(i, s, addTraj);
        }
        workQueue.enqueue_bulk(buffer, bCounter);

        impactWorker(assigned - 1, sPtr, addTraj);

        while(threadCount < assigned){
            std::this_thread::yield();
        }
        
        for(int i=0; i<assigned-1; i++){
            threads[i].join();
        }

        s.completedImpact = true;
    }

    //Post-Penetration Section
    private:
    double dtf = 0.0001;
    double xf0 = 0, yf0 = 0;

    void postPenTraj(const unsigned int i, shell& s, double v_x, double v_y, double v_z, double thickness){
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();
        double T, p, rho, t; 

        /* [indices]           0     1     2     3 [Bits]
         * pos                 x     y     z     - 256
         * velocities          v_x   v_y   v_z   - 256
         * velocitiesSquared   v_x^2 v_y^2 v_z^2 - 256
         * dragIntermediary    ad_x  ad_y  ad_z  - 256
         * xz_dragIntermediary ad_x  ad_z  -     - 128
         */
        __m256d pos, velocities, velocitiesSquared, dragIntermediary;
        __m128d xz_dragIntermediary; 
        pos[0] = xf0, pos[1] = yf0, pos[2] = xf0;
        velocities[0] = v_x, velocities[1] = v_y, velocities[2] = v_z;
        t = 0;
        if(v_x > 0){
            while(t < s.get_fuseTime()){
                pos = _mm256_fmadd_pd(velocities, _mm256_set1_pd(dtf), pos); //pos += velocities * dt
                
                //Calculate air density - likely unnecessary for this section as distances are so short
                T = t0 - L*pos[1];
                p = p0*pow((1-L*pos[1]/t0),(g*M/(R*L)));
                rho = p*M/(R*T);

                //Calculated drag deceleration
                velocitiesSquared = _mm256_mul_pd(velocities, velocities); //velocitiesSquared = velocities * velocities
                xz_dragIntermediary = _mm_mul_pd(_mm_set1_pd(k*rho),
                    _mm_fmadd_pd(_mm_set1_pd(cw_1), _mm_set_pd(velocitiesSquared[2], velocitiesSquared[0]), 
                    _mm_mul_pd(_mm_set1_pd(cw_2), _mm_set_pd(velocities[2], velocities[0]))
                    )
                ); //xz_dragIntermediary = (k * rho) * (cw_1 * velocitiesSquared[2, 0] + cw_2 * velocities[2, 0])
                
                dragIntermediary[0] = xz_dragIntermediary[0]; //x
                dragIntermediary[1] = (g - k*rho*(cw_1*velocitiesSquared[1]+cw_2*fabs(velocities[1]))*signum(velocities[1])); 
                dragIntermediary[2] = xz_dragIntermediary[1]; //z

                velocities = _mm256_fmadd_pd(_mm256_set1_pd(dtf * -1), dragIntermediary, velocities); //velocities -= dtf * dragIntermediary
                t += dtf;                                                                             
            }
            s.getPostPen(i, post::x) = pos[0];
            s.getPostPen(i, post::y) = pos[1];
            s.getPostPen(i, post::z) = pos[2];
            s.getPostPen(i, post::xwf) = (thickness >= s.get_threshold()) * pos[0] + !(thickness >= s.get_threshold()) * -1;
        }else{
            s.getPostPen(i, post::x) = 0;
            s.getPostPen(i, post::y) = 0;
            s.getPostPen(i, post::z) = 0;
            s.getPostPen(i, post::xwf) = 0;
        }
    }

    //template<typename T>
    void fillCopy(int id, shell* s, std::vector<double>* angles){
        //std::cout<<id<<"\n";
        for(int i=angles->size() * id / assigned; i<angles->size() * (id + 1) / assigned; i++){
            std::fill_n(s->postPenData.begin() + i * s->impactSize, s->impactSize, (double) angles->at(i));
            std::copy_n(s->getImpactPtr(0, impact::distance), s->impactSize, s->postPenData.begin() + s->postPenSize + i * s->impactSize);
        }
        counter.fetch_add(1, std::memory_order_relaxed);
    }

    //template<typename T>
    void parallelFillCopy(shell* s, std::vector<double>* angles){
        std::vector<std::thread> threads;
        counter = 0;
        length = angles->size();
        if(length < std::thread::hardware_concurrency()){
            assigned = length;
        }else{
            assigned = std::thread::hardware_concurrency();
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

    void calculatePostPen(const double thickness, shell& s, std::vector<double>& angles){

        if(!s.completedImpact){
            std::cout<<"Standard Not Calculated - Running automatically\n";
            calculateImpact(s, false);
        }

        s.postPenSize = s.impactSize * angles.size();
        s.postPenData.resize(6 * s.postPenSize);

        parallelFillCopy(&s, &angles);        

        std::vector<std::thread> threads;
        counter = 0, threadCount = 0;
        length = s.postPenSize/vSize;
        if(length < std::thread::hardware_concurrency()){
            assigned = length;
        }else{
            assigned = std::thread::hardware_concurrency();
        }
        int buffer[16];
        int bCounter = 0;
        for(int i=0; i<s.postPenSize; i+=vSize){
            buffer[bCounter] = i;
            bCounter++;
            if(bCounter == 16){
                workQueue.enqueue_bulk(buffer, 16);
                bCounter = 0;
            }
        }
        workQueue.enqueue_bulk(buffer, bCounter);

        shell* ptr = &s;
        for(int i=0; i<assigned - 1; i++){
            threads.emplace_back([=]{postPenWorker(i, thickness, ptr);});
        }
        postPenWorker(assigned - 1, thickness, &s);

        while(threadCount < assigned){
            std::this_thread::yield();
        }

        for(int i=0; i<assigned - 1; i++){
            threads[i].join();
        }
        s.completedPostPen = true;
        
    }

    void postPenWorker(int thread, double thickness, shell* s){
        while(counter < length){
            int index;
            if(workQueue.try_dequeue(index)){
                multiPostPen(index, thickness, *s);
                counter.fetch_add(1, std::memory_order_relaxed);
            }
            else{
                std::this_thread::yield();
            }
        }
        threadCount.fetch_add(1, std::memory_order_relaxed);
    }
    
    void multiPostPen(int i, const double thickness, shell& s){
        //std::cout<<index<<"\n";
        //unsigned int i = index * vSize;
        __m256d hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV;
        __m256d v0V, pPVV, ePenetrationV, eThickness, hFAngleV, vFAngleV;
        __m256d v_x, v_y, v_z;
        unsigned int distIndex = (i < s.impactSize) ? i : i % s.impactSize;
        unsigned int anglesIndex = i / s.impactSize;

        unsigned int j, k = 0;

        if(i + vSize <= s.postPenSize){
            hAngleV = _mm256_loadu_pd(&s.postPenData[i]);
        }else{
            for(j = 0; (i + j)< s.postPenSize; j++){
                hAngleV[j] = s.postPenData[i+j];
            }
        }

        if(distIndex < s.impactSize - vSize + 1){
            vAngleV = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::impactAHR));
            ePenetrationV = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::rawPen));
            v0V = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::impactV));
        }else{
            for(j = 0; (j + distIndex < s.impactSize) && (j < vSize); j++){
                vAngleV[j] = s.getImpact(distIndex + j, impact::impactAHR);
                ePenetrationV[j] = s.getImpact(distIndex + j, impact::rawPen);
                v0V[j] = s.getImpact(distIndex + j, impact::impactV);
            }
            if(anglesIndex < s.postPenSize / s.impactSize){
                for(; (j < vSize); j++){
                    vAngleV[j] = s.getImpact(k, impact::impactAHR);
                    ePenetrationV[j] = s.getImpact(k, impact::rawPen);
                    v0V[j] = s.getImpact(k, impact::impactV);
                    k++;
                }
            }
        }

        hAngleV = _mm256_mul_pd(hAngleV, _mm256_set1_pd(M_PI/180));            
        cAngleV = xacos(_mm256_mul_pd(xcos(hAngleV), xcos(vAngleV)));
        nCAngleV = calcNormalizationRSIMD(cAngleV, s.get_normalizationR());
        eThickness = _mm256_div_pd(_mm256_set1_pd(thickness), xcos(nCAngleV));
            
        pPVV = _mm256_max_pd(_mm256_mul_pd(v0V, 
            _mm256_sub_pd(_mm256_set1_pd(1), 
                xexp(_mm256_sub_pd(_mm256_set1_pd(1),
                    _mm256_div_pd(ePenetrationV, eThickness)
                )))
        ), _mm256_set1_pd(0));

        aAngleV = _mm256_div_pd(nCAngleV, cAngleV);
        hFAngleV = _mm256_mul_pd(hAngleV, aAngleV);
        vFAngleV = _mm256_mul_pd(vAngleV, aAngleV);
        
        __m256d vFAngleVCos = xcos(vFAngleV);
        v_x = _mm256_mul_pd(pPVV,
            _mm256_mul_pd(vFAngleVCos, xcos(hFAngleV))
        );
        v_y = _mm256_mul_pd(pPVV, xsin(vFAngleV));
        v_z = _mm256_mul_pd(pPVV,
            _mm256_mul_pd(vFAngleVCos, xsin(hFAngleV))
        );

        for(unsigned int j=0; (j<vSize) && (j+i < s.postPenSize); j++){
            postPenTraj(i+j, s, v_x[j], v_y[j], v_z[j], eThickness[j]);
        }
        //std::cout<<index<<" Completed\n";
    }
};

}