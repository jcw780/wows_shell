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

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <cstring>
#include <algorithm>
#include <type_traits>
		
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
enum stdDataIndex{
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
    unsigned int size, postPenSize; //number of distances in: standard, postPen
    unsigned int sizeAligned, postPenSizeAligned; //Not 100% necessary - sizes adjusted to fulfill alignment
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
    std::vector<double> postPenData;
    
    double& getImpact(unsigned int i, unsigned int j){
        return impactData[i + j * sizeAligned];
    }

    double* getImpactPtr(unsigned int i, unsigned int j){
        return impactData.data() + i + j * sizeAligned;
    }

    double& getPostPen(unsigned int i, unsigned int j){
        return postPenData[i + j * postPenSize];
    }

    double* getPostPenPtr(unsigned int i, unsigned int j){
        return postPenData.data() + i + j * postPenSize;
    }


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
            for(int j=0; j<post::maxColumns; j++){
                std::cout<< std::fixed<< std::setprecision(4) << postPenData[i + j * postPenSize] << " "; 
            }
            std::cout<<"\n";
        }
        std::cout<<"Completed Post-Penetration\n";
    }

    void printImpactData(){
        for(unsigned int i=0; i<size; i++){
            for(unsigned int j=0; j<impact::maxColumns; j++){
                std::cout<<std::fixed<<std::setprecision(4)<<getImpact(i, j)<< " ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"Completed Standard Data"<<std::endl;
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
    //Physical Constants
    //                    Description                    Units
    double g = 9.81;      //Gravitational Constant       m/(s^2)
    double t0 = 288;      //Temperature at Sea Level     K
    double L = 0.0065;    //Atmospheric Lapse Rate       C/m
    double p0 = 101325;   //Pressure at Sea Level        Pa
    double R = 8.31447;   //Ideal Gas Constant           J/(mol K)
    double M = 0.0289644; //Molarity of Air at Sea Level kg/mol
    double cw_1 = 1;
    
    //Calculation Parameters Description          Units
    double max = 25;         //Max Angle          degrees
    double min = 0;          //Min Angle          degrees
    double precision = .1;   //Angle Step         degrees
    double x0 = 0, y0 = 0;   //Starting x0, y0    m
    double dt = .01;         //Time step          s

    static_assert(sizeof(double) == 8, "Size of double is not 8 - required for AVX2"); //Use float64 in the future
    static_assert(std::numeric_limits<double>::is_iec559, "Type is not IEE754 compliant");
    static constexpr unsigned int vSize = (256 / 8) / sizeof(double);
    
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
        //s.stdData[i+j+s.sizeAligned*impact::distance] = pos[0];
        s.getImpact(i + j, impact::distance) = pos[0];
        vx[j] = velocity[0];
        vy[j] = velocity[1];
        tSIMD[j] = t;
    }

    void quadTraj(const unsigned int i, shell& s, const bool addTraj){
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

        for(unsigned int j = 0; (j+i<s.size) && (j < vSize); j++){
            singleTraj(i, j, s, vx, vy, tSIMD);
        }

        //Calculate [2]Impact Angles (impactAHR) , [7] Impact Angle Deck (impactADD)
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

        //Calculate [3]Impact Velocity (impactV),  [4]Raw Penetration (rawPen)
        iVSIMD = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(vx, vx), _mm256_mul_pd(vy, vy)));
        _mm256_storeu_pd(s.getImpactPtr(i, impact::impactV), iVSIMD);

        rPSIMD = _mm256_mul_pd(xpow(iVSIMD, _mm256_set1_pd(1.1)), _mm256_set1_pd(pPPC));
        _mm256_storeu_pd(s.getImpactPtr(i, impact::rawPen), rPSIMD);

        //Calculate [5]EPH  [8]EPV
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenH), _mm256_mul_pd(xcos(angleRSIMD),rPSIMD));

        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenD), _mm256_mul_pd(xcos(angleSIMD),rPSIMD));
        
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenHN),
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleRSIMD, normalizationR))
            )
        );
        _mm256_storeu_pd(s.getImpactPtr(i, impact::ePenDN),
            _mm256_mul_pd(
                rPSIMD, xcos(calcNormalizationRSIMD(angleSIMD, normalizationR))
            )
        );

        _mm256_storeu_pd(s.getImpactPtr(i, impact::tToTarget), tSIMD);

        _mm256_storeu_pd(s.getImpactPtr(i, impact::tToTargetA), _mm256_div_pd(tSIMD, 
            _mm256_set1_pd(3.1))
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
        unsigned int i;
        s.size = (unsigned int) (max - min) / precision;
        s.sizeAligned = (sizeof(__m256d)/sizeof(double) - (s.size % (sizeof(__m256d)/sizeof(double)))) + s.size;
        //s.sizeAligned = s.size;

        s.trajectories.resize(2 * s.size);
        s.impactData.resize(impact::maxColumns * s.sizeAligned);

        //omp_set_num_threads(6);
        #pragma omp parallel for schedule(dynamic, 2)
        for(i=0; i<s.size; i+=vSize){
            quadTraj(i, s, addTraj);
        }
        s.completedImpact = true;
    }

    //Post-Penetration Section
    private:
    double dtf = 0.0001;
    double xf0 = 0, yf0 = 0;
    //bool completed = false;

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
            //std::cout<<i<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"\n";
            s.postPenData[i+s.postPenSize*2] = pos[0];
            s.postPenData[i+s.postPenSize*3] = pos[1];
            s.postPenData[i+s.postPenSize*4] = pos[2];
            if(thickness >= s.get_threshold()){
                s.postPenData[i+s.postPenSize*5] = pos[0];
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

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE SIMD ALIGNMENT MAY NOT BE GUARANTEED
     * [0:1) Lateral Angle [1:2) Distance [2:3) X [3:4) Y [4:5) Z [5:6) XWF
     */
    //std::vector<double, AlignmentAllocator<double, 32>> postPenData;

    template<typename T>
    void calculatePostPen(const double thickness, shell& s, std::vector<T>& angles){
        static_assert(std::is_arithmetic<T>(), "Cannot use non numeric type");

        if(!s.completedImpact){
            std::cout<<"Standard Not Calculated - Running automatically\n";
            calculateImpact(s, false);
        }

        s.postPenSize = s.size * angles.size();
        s.postPenData.resize(6 * s.postPenSize);


        #pragma omp parallel for
        for(unsigned int i=0; i < angles.size(); i++){
            std::fill_n(s.postPenData.begin() + i * s.size, s.size, (double) angles[i]);
            std::copy_n(s.getImpactPtr(0, impact::distance), s.size, s.postPenData.begin() + s.postPenSize + i * s.size);
        }

        #pragma omp parallel for
        for(unsigned int i=0; i < s.postPenSize; i+=vSize){
            __m256d hAngleV, vAngleV, cAngleV, nCAngleV, aAngleV;
            __m256d v0V, pPVV, ePenetrationV, eThickness, hFAngleV, vFAngleV;
            __m256d v_x, v_y, v_z;
            unsigned int distIndex = (i < s.size) ? i : i % s.size;
            unsigned int anglesIndex = i / s.size;

            unsigned int j, k = 0;

            if(i + vSize <= s.postPenSize){
                hAngleV = _mm256_loadu_pd(&s.postPenData[i]);
            }else{
                for(j = 0; (i + j)< s.postPenSize; j++){
                    hAngleV[j] = s.postPenData[i+j];
                }
            }

            if(distIndex < s.size - vSize + 1){
                vAngleV = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::impactAHR));
                ePenetrationV = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::rawPen));
                v0V = _mm256_loadu_pd(s.getImpactPtr(distIndex, impact::impactV));
            }else{
                for(j = 0; (j + distIndex < s.size) && (j < vSize); j++){
                    vAngleV[j] = s.getImpact(distIndex + j, impact::impactAHR);
                    ePenetrationV[j] = s.getImpact(distIndex + j, impact::rawPen);
                    v0V[j] = s.getImpact(distIndex + j, impact::impactV);
                }
                if(anglesIndex < angles.size()){
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
        }
        s.completedPostPen = true; 
    }


};

}