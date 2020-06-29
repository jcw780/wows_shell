#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

//#include "concurrentqueue/concurrentqueue.h"
#include "utility.hpp"


namespace shell {

namespace impact {
static constexpr unsigned int maxColumns = 13;
static constexpr unsigned int maxColumnsFit = 7;
enum impactDataIndex {
    distance,
    launchAngle,
    impactAngleHorizontalRadians,  // Negative for Falling
    impactAngleHorizontalDegrees,  // Positive for Falling
    impactVelocity,
    timeToTarget,
    timeToTargetAdjusted,  // Adjusted for in game shell time
    rawPenetration,
    effectivePenetrationHorizontal,
    effectivePenetrationHorizontalNormalized,
    impactAngleDeckDegrees,
    effectivePenetrationDeck,
    effectivePenetrationDeckNormalized,
};
static_assert(effectivePenetrationDeckNormalized == (maxColumns - 1),
              "Invalid standard columns");
}  // namespace impact

namespace angle {
static constexpr unsigned int maxColumns = 9;
enum angleDataIndex {
    distance,
    ricochetAngle0Radians,
    ricochetAngle0Degrees,
    ricochetAngle1Radians,
    ricochetAngle1Degrees,
    armorRadians,
    armorDegrees,
    fuseRadians,
    fuseDegrees
};
static_assert(fuseDegrees == (maxColumns - 1), "Invalid angle columns");
}  // namespace angle

namespace post {
static constexpr unsigned int maxColumns = 6;
enum postPenDataIndex { angle, distance, x, y, z, xwf };
static_assert(xwf == (maxColumns - 1), "Invaild postpen columns");
}  // namespace post

enum numerical { forwardEuler, rungeKutta2, rungeKutta4, adamsBashforth5 };

class shell {
   public:                 // Description                units
    double v0;             // muzzle velocity            m/s
    double caliber;        // shell caliber              m
    double krupp;          // shell krupp                [ndim]
    double mass;           // shell mass                 kg
    double cD;             // coefficient of drag        [ndim]
    double normalization;  // shell normalization        degrees
    double threshold;      // fusing threshold           mm
    double fuseTime;       // fuse time                  s
    double ricochet0;      // ricochet angle 0           degrees
    double ricochet1;      // ricochet angle 1           degrees
    double nonAP;          // penetration of non ap ammo mm
    bool enableNonAP;
    std::string name;

    double k, cw_2, pPPC, normalizationR, ricochet0R, ricochet1R;

    // Condenses initial values into values used by calculations
    //[Reduces repeated computations]
   public:
    void preProcess() {
        k = 0.5 * cD * pow((caliber / 2), 2) * M_PI / mass;
        // condensed drag coefficient
        cw_2 = 100 + 1000 / 3 * caliber;
        // quadratic drag coefficient
        pPPC = 0.5561613 * krupp / 2400 * pow(mass, 0.5506) /
               pow((caliber * 1000), 0.6521);
        // condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI;
        // normalization (radians)
        ricochet0R = ricochet0 / 180 * M_PI;
        ricochet1R = ricochet1 / 180 * M_PI;
    }

    unsigned int impactSize,
        postPenSize;  // number of distances in: standard, postPen
    unsigned int impactSizeAligned, postPenSizeAligned;
    // Not 100% necessary - sizes adjusted to fulfill alignment
    bool completedImpact = false, completedAngles = false,
         completedPostPen = false;

    /*trajectories output
    [0           ]trajx 0        [1           ]trajy 1
    ...
    [size * 2 - 2]trajx size - 1 [size * 2 - 1]trajy size - 1
    */
    std::vector<std::vector<double>> trajectories;

    // Refer to stdDataIndex enums defined above
    std::vector<double> impactData;

    /* Angles data
     * [0:1)-ra0 max lateral angle
     * [1:2)-ra0 max lateral angle
     * [2:3)-ra1 max lateral angle
     * [3:4)-ra1 max lateral angle
     * [4:5)-penetration max lateral angle
     * [5:6)-penetration max lateral angle
     */
    std::vector<double> angleData;

    /* WARNING: LOCATION OF LATERAL ANGLE IN VECTOR CANNOT BE CHANGED OR ELSE
     * SIMD ALIGNMENT MAY NOT BE GUARANTEED [0:1) Lateral Angle [1:2) Distance
     * [2:3) X [3:4) Y [4:5) Z             [5:6) XWF See enums defined above
     */
    std::vector<double> postPenData;

    shell() = default;

    shell(const double caliber, const double v0, const double cD,
          const double mass, const double krupp, const double normalization,
          const double fuseTime, const double threshold, const double ricochet0,
          const double ricochet1, const double nonAP, const std::string &name) {
        setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                  threshold, ricochet0, ricochet1, nonAP, name);
    }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const double nonAP,
                   const std::string &name) {
        //                                                    Impact   PostPen
        this->fuseTime = fuseTime;  // Shell fusetime        | No     | Yes
        this->v0 = v0;              // Shell muzzle velocity | Yes    | No
        this->caliber = caliber;    // Shell caliber         | Yes    | Yes
        this->krupp = krupp;        // Shell krupp           | Yes    | Yes
        this->mass = mass;          // Shell mass            | Yes    | Yes
        this->normalization =
            normalization;  // Shell normalization   | Yes    | Yes
        this->cD = cD;      // Shell air drag coefficient | Yes    | Yes
        this->name = name;  // Shell name  | No     | No
        this->ricochet0 = ricochet0;  // Ricochet Angle 0 | No     | Yes
        this->ricochet1 = ricochet1;  // Ricochet Angle 1 | No     | Yes
                                      // Shell fusing threshold | No     | Yes
        this->threshold = threshold;
        this->nonAP = nonAP;
        if (this->nonAP > 0) {
            this->enableNonAP = true;
        } else {
            this->enableNonAP = false;
        }
        preProcess();
    }

    // Getter Functions
    double &get_impact(unsigned int row, unsigned int impact) {
        return impactData[row + impact * impactSizeAligned];
    }

    double *get_impactPtr(unsigned int row, unsigned int impact) {
        return impactData.data() + row + impact * impactSizeAligned;
    }

    double &get_angle(unsigned int row, unsigned int impact) {
        return angleData[row + impact * impactSizeAligned];
    }

    double *get_anglePtr(unsigned int row, unsigned int impact) {
        return angleData.data() + row + impact * impactSizeAligned;
    }

    double &get_postPen(unsigned int row, unsigned int angle,
                        unsigned int impact) {
        return postPenData[row + angle * postPenSize + impact * impactSize];
    }

    double *get_postPenPtr(unsigned int row, unsigned int angle,
                           unsigned int impact) {
        return postPenData.data() + row + angle * postPenSize +
               impact * impactSize;
    }

    // Linear interpolate by distance - WIP

    double interpolateDistanceImpact(double distance, unsigned int impact) {
        auto iter_max = std::lower_bound(
            get_impactPtr(0, impact::impactDataIndex::distance),
            get_impactPtr(0, impact::impactDataIndex::distance + 1), distance);
        double maxDistance = *iter_max;
        unsigned int maxIndex =
            iter_max - get_impactPtr(0, impact::impactDataIndex::distance);
        double maxTarget = get_impact(maxIndex, impact);

        auto iter_min = iter_max - 1;
        double minDistance = *iter_min;
        unsigned int minIndex =
            iter_min - get_impactPtr(0, impact::impactDataIndex::distance);
        double minTarget = get_impact(minIndex, impact);

        /*std::cout << minIndex << " " << minDistance << " " << minTarget << " "
                  << maxIndex << " " << maxDistance << " " << maxTarget << "\n";
                  */
        double slope = ((maxTarget - minTarget) / (maxDistance - minDistance));
        return slope * (distance - minDistance) + minTarget;
    }

    // internal computed data - fixed
    double get_v0() { return v0; }
    double get_k() { return k; }
    double get_cw_2() { return cw_2; }
    double get_pPPC() { return pPPC; }
    double get_normalizationR() { return normalizationR; }

    void printAngleData() {
        for (unsigned int i = 0; i < impactSize; i++) {
            for (unsigned int j = 0; j < angle::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_angle(i, j) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Angle Data\n";
    }

    void printPostPenData() {
        for (unsigned int i = 0; i < postPenSize; i++) {
            for (unsigned int j = 0; j < post::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_postPen(i, j, 0) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Completed Post-Penetration\n";
    }

    void printImpactData() {
        for (unsigned int i = 0; i < impactSize; i++) {
            for (unsigned int j = 0; j < impact::maxColumns; j++) {
                std::cout << std::fixed << std::setprecision(4)
                          << get_impact(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Completed Standard Data" << std::endl;
    }
    void printTrajectory(unsigned int target) {
        if (target >= impactSize) {
            std::cout << "Target Not Within Range of: " << impactSize
                      << std::endl;
        } else {
            printf("Index:[%d] X Y\n", target);
            for (std::vector<double>::size_type i = 0;
                 i < trajectories[target * 2].size(); i++) {
                std::cout << trajectories[target * 2][i] << " "
                          << trajectories[target * 2 + 1][i] << std::endl;
            }
        }
    }
};

class shellCalc {
   private:
    // Physical Constants       Description                  | Units
    double g = 9.81;       // Gravitational Constant       | m/(s^2)
    double t0 = 288;       // Temperature at Sea Level     | K
    double L = 0.0065;     // Atmospheric Lapse Rate       | C/m
    double p0 = 101325;    // Pressure at Sea Level        | Pa
    double R = 8.31447;    // Ideal Gas Constant           | J/(mol K)
    double M = 0.0289644;  // Molarity of Air at Sea Level | kg/mol
    double cw_1 = 1;

    double gMRL = (g * M / (R * L));
    // Calculation Parameters
    double max = 25;        // Max Angle                    | degrees
    double min = 0;         // Min Angle                    | degrees
    double precision = .1;  // Angle Step                   | degrees
    double x0 = 0, y0 = 0;  // Starting x0, y0              | m
    double dt_min = .02;    // Time step                    | s

    // delta t (dtf) for fusing needs to be smaller than the delta t (dt) used
    // for trajectories due to the shorter distances. Otherwise results become
    // jagged - precision suffers.
    double dtf = 0.0001;
    double xf0 = 0, yf0 = 0;

    // For vectorization - though probably not 100% necessary anymore since
    // intrinsics were removed [intrinsics had no significant improvements in
    // runtime] - but might still make it easier to vectorize
    static_assert(
        sizeof(double) == 8,
        "Size of double is not 8 - required for AVX2");  // Use float64
                                                         // in the future
    static_assert(std::numeric_limits<double>::is_iec559,
                  "Type is not IEE754 compliant");
    static constexpr unsigned int vSize = (256 / 8) / sizeof(double);
    static constexpr unsigned int minTasksPerThread = vSize;

   public:
    double calcNormalizationR(
        const double angle,
        const double normalizationR) {  // Input in radians
        return (fabs(angle) > normalizationR) * (fabs(angle) - normalizationR);
    }

    inline int signum(double x) { return ((0.0) < x) - (x < (0.0)); }

    shellCalc() = default;
    // Replace with setter in the future
    void editTestParameters(const double max, const double min,
                            const double precision, const double x0,
                            const double y0, const double dt_min,
                            const double xf0, const double yf0,
                            const double dtf) {
        this->max = max;
        this->min = min;
        this->precision = precision;
        this->x0 = x0;
        this->y0 = y0;
        this->dt_min = dt_min;
        this->xf0 = xf0;
        this->yf0 = yf0;
        this->dtf = dtf;
    }

    void set_max(const double max) { this->max = max; }
    void set_min(const double min) { this->min = min; }
    void set_precision(const double precision) { this->precision = precision; }
    void set_x0(const double x0) { this->x0 = x0; }
    void set_y0(const double y0) { this->y0 = y0; }
    void set_dt_min(const double dt) { this->dt_min = dt; }
    void set_xf0(const double xf0) { this->xf0 = xf0; }
    void set_yf0(const double yf0) { this->yf0 = yf0; }
    void set_dtf(const double dtf) { this->dtf = dtf; }

   private:
    // Utility functions
    // mini 'threadpool' used to kick off multithreaded functions
    template <typename...>
    inline static constexpr bool dependent_false = false;

    template <typename O, typename F, typename... Args>
    void mtFunctionRunner(int assigned, int length, int size, O object,
                          F function, Args... args) {
        if (assigned > 1) {
            mtFunctionRunnerSelected<true>(assigned, length, size, object,
                                           function, args...);
        } else {
            mtFunctionRunnerSelected<false>(assigned, length, size, object,
                                            function, args...);
        }
    }

    static constexpr int workQueueBufferSize = 16;
    template <bool multiThreaded, typename O, typename F, typename... Args>
    void mtFunctionRunnerSelected(unsigned int assigned, unsigned int length,
                                  unsigned int size, O object, F function,
                                  Args... args) {
        if constexpr (multiThreaded) {
            std::atomic<int> counter{0};
            //moodycamel::ConcurrentQueue<int> workQueue;
            std::vector<std::thread> threads(assigned - 1);
            for (unsigned int i = 0; i < assigned - 1; i++) {
                threads[i] = std::thread([&] {
                    mtWorker(counter, length, i, object, function,
                             args...);
                });
            }

            mtWorker(counter, length, assigned - 1, object, function,
                     args...);

            for (unsigned int i = 0; i < assigned - 1; i++) {
                threads[i].join();
            }
        } else {
            for (unsigned int i = 0; i < size; i += vSize) {
                (object->*function)(i, args...);
            }
        }
    }

    template <typename O, typename F, typename... Args>
    void mtWorker(std::atomic<int> &counter, const int length,
                  const int threadID, O object, F function, Args... args) {
        // threadID is largely there for debugging
        while (counter < length) {
            int index = counter.fetch_add(1, std::memory_order_relaxed);
            if(index < length){
                //std::cout<<index<<"\n";
                (object->*function)(index * vSize, args...);
            }
        }
    }

    unsigned int assignThreadNum(unsigned int length, unsigned int nThreads) {
        if (length > nThreads * minTasksPerThread) {
            return nThreads;
        } else {
            return ceil((double)length / minTasksPerThread);
        }
    }

    // Numerical Analysis Methods and Dependencies
    // Calculate changes given conditions
    template <unsigned int dims>
    void calcDeltas(const double *const current, double *const deltas,
                    const double dt, const double k, const double cw_2) {
        enum position { x, y, z };
        enum velocity { v_x, v_y, v_z };

        // std::array<double, dims * 2> deltas;
        for (unsigned int i = 0; i < dims; i++) {
            deltas[i] = current[i + dims] * dt;
            // currentL[i] += deltas[i];
        }

        /*for (unsigned int i=0; i<dims * 2; i++){
            std::cout<<current[i]<<" ";
        }*/

        double T, p, rho;
        double postChangeY = current[position::y] + deltas[position::y];
        // Calculating air density
        T = t0 - L * postChangeY;
        // Calculating air temperature at altitude
        p = p0 * pow((1 - L * (postChangeY) / t0), gMRL);
        // Calculating air pressure at altitude
        rho = p * M / (R * T);
        // Use ideal gas law to calculate air density
        // std::cout<<" rho "<<rho<<" ";
        // Calculate drag deceleration
        std::array<double, dims> velocitySquared;
        for (unsigned int l = 0; l < dims; l++) {
            velocitySquared[l] = current[dims + l] * current[dims + l];
        }  // v^2 = v * v

        for (unsigned int i = 0; i < dims; i++) {
            deltas[dims + i] = cw_1 * velocitySquared[i];
        }

        deltas[dims + velocity::v_x] += cw_2 * current[dims + velocity::v_x];
        deltas[dims + velocity::v_y] +=
            cw_2 * fabs(current[dims + velocity::v_y]);
        deltas[dims + velocity::v_y] *= signum(current[dims + velocity::v_y]);

        if constexpr (dims >= 3) {
            deltas[dims + velocity::v_z] +=
                cw_2 * current[dims + velocity::v_z];
        }

        for (unsigned int i = 0; i < dims; i++) {
            deltas[dims + i] *= (k * rho);
        }

        deltas[dims + velocity::v_y] = g + deltas[dims + velocity::v_y];
        for (unsigned int i = 0; i < dims; i++) {  // v -= drag * dt
            deltas[dims + i] *= (-1 * dt);
        }
    }

    // Numerical Analysis Methods
    template <unsigned int dims>
    void forwardEuler(std::array<double, 2 * dims> &input, const double dt,
                      const double k, const double cw_2) {
        constexpr unsigned int arrSize = 2 * dims;
        std::array<double, arrSize> delta;
        calcDeltas<dims>(input.data(), delta.data(), dt, k, cw_2);
        for (unsigned int i = 0; i < arrSize; i++) {
            input[i] += delta[i];
        }
    }

    template <unsigned int dims>
    void rungeKutta2(std::array<double, 2 * dims> &input, const double dt,
                     const double k, const double cw_2) {
        constexpr unsigned int arrSize = 2 * dims;
        std::array<double, 3 * arrSize> flatArr;
        enum fAI { k1I, k2I, intermediateI };
        static_assert((fAI::intermediateI + 1) * arrSize == flatArr.size());
        double *k1 = &flatArr[fAI::k1I * arrSize];
        double *k2 = &flatArr[fAI::k2I * arrSize];
        double *intermediate = &flatArr[fAI::intermediateI * arrSize];

        calcDeltas<dims>(input.data(), k1, dt, k, cw_2);
        utility::fmaArr<arrSize>(.5, k1, input.data(), intermediate);
        calcDeltas<dims>(intermediate, k2, dt, k, cw_2);
        utility::addArrInplace<arrSize>(k2, input.data());
    }

    template <unsigned int dims>
    void rungeKutta4(std::array<double, 2 * dims> &input, const double dt,
                     const double k, const double cw_2) {
        constexpr unsigned int arrSize = 2 * dims;
        std::array<double, 5 * arrSize> flatArr;
        enum fAI { k1I, k2I, k3I, k4I, intermediateI };
        static_assert((intermediateI + 1) * arrSize == flatArr.size());
        double *k1 = &flatArr[fAI::k1I * arrSize];
        double *k2 = &flatArr[fAI::k2I * arrSize];
        double *k3 = &flatArr[fAI::k3I * arrSize];
        double *k4 = &flatArr[fAI::k4I * arrSize];
        double *intermediate = &flatArr[fAI::intermediateI * arrSize];

        calcDeltas<dims>(input.data(), k1, dt, k, cw_2);
        utility::fmaArr<arrSize>(.5, k1, input.data(), intermediate);
        calcDeltas<dims>(intermediate, k2, dt, k, cw_2);

        utility::fmaArr<arrSize>(.5, k2, input.data(), intermediate);
        calcDeltas<dims>(intermediate, k3, dt, k, cw_2);

        utility::fmaArr<arrSize>(1, k3, input.data(), intermediate);
        calcDeltas<dims>(intermediate, k4, dt, k, cw_2);

        utility::fmaArrInplace<arrSize>(2, k2, k1);
        utility::fmaArrInplace<arrSize>(2, k3, k1);
        utility::fmaArrInplace<arrSize>(1, k4, k1);
        utility::fmaArrInplace<arrSize>((1.0 / 6.0), k1, input.data());
    }

    template <bool addTraj, unsigned int dims, typename Comparision>
    unsigned short adamsBashforth5(Comparision comp, std::vector<double> &tx,
                                   std::vector<double> &ty,
                                   std::array<double, 2 * dims> &input,
                                   const double dt, const double k,
                                   const double cw_2) {
        constexpr unsigned int arrSize = 2 * dims;
        std::array<double, arrSize * 5> flatArr;
        enum fAI { d0I, d1I, d2I, d3I, d4I };
        static_assert((d4I + 1) * arrSize == flatArr.size());
        double *d0 = &flatArr[d0I * arrSize];
        double *d1 = &flatArr[d1I * arrSize];
        double *d2 = &flatArr[d2I * arrSize];
        double *d3 = &flatArr[d3I * arrSize];
        double *d4 = &flatArr[d4I * arrSize];

        unsigned short counter = 0;
        calcDeltas<dims>(input.data(), d0, dt, k, cw_2);
        utility::addArrInplace<arrSize>(d0, input.data());
        counter++;
        if constexpr (addTraj) {
            tx.push_back(input[singleTrajVars::x]);
            ty.push_back(input[singleTrajVars::y]);
        }

        if (!(this->*comp)(input, dt)) {
            return counter;
        }

        calcDeltas<dims>(input.data(), d1, dt, k, cw_2);
        utility::fmaArrInplace<arrSize>((3.0 / 2.0), d1, input.data());
        utility::fmaArrInplace<arrSize>((-1.0 / 2.0), d0, input.data());
        counter++;
        if constexpr (addTraj) {
            tx.push_back(input[singleTrajVars::x]);
            ty.push_back(input[singleTrajVars::y]);
        }

        if (!(this->*comp)(input, dt)) {
            return counter;
        }

        calcDeltas<dims>(input.data(), d2, dt, k, cw_2);
        utility::fmaArrInplace<arrSize>((23.0 / 12.0), d2, input.data());
        utility::fmaArrInplace<arrSize>((-16.0 / 12.0), d1, input.data());
        utility::fmaArrInplace<arrSize>((5.0 / 12.0), d0, input.data());
        counter++;
        if constexpr (addTraj) {
            tx.push_back(input[singleTrajVars::x]);
            ty.push_back(input[singleTrajVars::y]);
        }

        if (!(this->*comp)(input, dt)) {
            return counter;
        }

        calcDeltas<dims>(input.data(), d3, dt, k, cw_2);
        utility::fmaArrInplace<arrSize>((55.0 / 24.0), d3, input.data());
        utility::fmaArrInplace<arrSize>((-59.0 / 24.0), d2, input.data());
        utility::fmaArrInplace<arrSize>((37.0 / 24.0), d1, input.data());
        utility::fmaArrInplace<arrSize>((-9.0 / 24.0), d0, input.data());
        counter++;
        if constexpr (addTraj) {
            tx.push_back(input[singleTrajVars::x]);
            ty.push_back(input[singleTrajVars::y]);
        }

        if (!(this->*comp)(input, dt)) {
            return counter;
        }

        calcDeltas<dims>(input.data(), d4, dt, k, cw_2);
        utility::fmaArrInplace<arrSize>((1901.0 / 720.0), d4, input.data());
        utility::fmaArrInplace<arrSize>((-2774.0 / 720.0), d3, input.data());
        utility::fmaArrInplace<arrSize>((2616.0 / 720.0), d2, input.data());
        utility::fmaArrInplace<arrSize>((-1274.0 / 720.0), d1, input.data());
        utility::fmaArrInplace<arrSize>((251.0 / 720.0), d0, input.data());
        counter++;
        if constexpr (addTraj) {
            tx.push_back(input[singleTrajVars::x]);
            ty.push_back(input[singleTrajVars::y]);
        }

        return counter;
    }

    // Please edit numerical enum before changing this
    template <unsigned int Numerical, unsigned int Dims>
    void numericalMethod(std::array<double, 2 * Dims> &input, const double dt,
                         const double k, const double cw_2) {
        if constexpr (Numerical == numerical::forwardEuler) {
            forwardEuler<Dims>(input, dt, k, cw_2);
        } else if (Numerical == numerical::rungeKutta2) {
            rungeKutta2<Dims>(input, dt, k, cw_2);
        } else if (Numerical == numerical::rungeKutta4) {
            rungeKutta4<Dims>(input, dt, k, cw_2);
        } else {
            // static_assert(dependent_false<Dims>, "Missing Singlestep
            // Function");
        }
    }

    template <bool addTraj, unsigned int Numerical, unsigned int Dims,
              typename Comparision>
    unsigned int multistepMethod(Comparision comp, std::vector<double> &tx,
                                 std::vector<double> &ty,
                                 std::array<double, 4> &input, const double dt,
                                 const double k, const double cw_2) {
        if constexpr (Numerical == numerical::adamsBashforth5) {
            return adamsBashforth5<addTraj, Dims>(comp, tx, ty, input, dt, k,
                                                  cw_2);
        } else {
            static_assert(dependent_false<Comparision>,
                          "Missing Multistep Function");
        }
    }

    template <unsigned int Numerical>
    static constexpr bool isMultistep() {
        if constexpr (Numerical == numerical::adamsBashforth5) {
            return true;
        } else {
            return false;
        }
    }

    // Impact Calculations Region
    enum singleTrajVars { x, y, v_x, v_y };
    static constexpr unsigned int singleTrajDims = 2;

    template <bool addTraj, unsigned int Numerical, unsigned int Dims,
              typename Comparision>
    void singleTrajStep(Comparision comp, std::vector<double> &tx,
                        std::vector<double> &ty, std::array<double, 4> &input,
                        double &t, const double dt, const double k,
                        const double cw_2) {
        if constexpr (isMultistep<Numerical>()) {
            t += (dt * multistepMethod<addTraj, Numerical, Dims>(
                           comp, tx, ty, input, dt, k, cw_2));
        } else {
            numericalMethod<Numerical, Dims>(input, dt, k, cw_2);
            t += dt;
            if constexpr (addTraj) {
                tx.push_back(input[singleTrajVars::x]);
                ty.push_back(input[singleTrajVars::y]);
            }
        }
    }

    bool stopFunction(std::array<double, singleTrajDims * 2> variables,
                      double dt) {
        return variables[singleTrajVars::y] >= 0;
    }

    bool stopFunctionHybrid(std::array<double, singleTrajDims * 2> variables,
                            double dt) {
        return variables[singleTrajVars::y] +
                   variables[singleTrajVars::v_y] * dt - g / 2 * dt * dt >=
               0;
    }

    template <bool AddTraj, unsigned int Numerical, bool Hybrid>
    void singleTraj(const unsigned int i, const unsigned int j, shell &s,
                    double *vx, double *vy, double *tVec) {
        static constexpr unsigned int __TrajBuffer__ = 128;
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();

        // unsigned int counter = 0;
        if constexpr (AddTraj) {
            s.trajectories[2 * (i + j)].clear();
            s.trajectories[2 * (i + j) + 1].clear();

            if (s.trajectories[2 * (i + j)].capacity() < __TrajBuffer__) {
                s.trajectories[2 * (i + j)].reserve(__TrajBuffer__);
            }

            if (s.trajectories[2 * (i + j) + 1].capacity() < __TrajBuffer__) {
                s.trajectories[2 * (i + j) + 1].reserve(__TrajBuffer__);
            }
        }
        // setting initial values
        std::array<double, 4> variables{x0, y0, vx[j], vy[j]};
        if constexpr (AddTraj) {
            s.trajectories[2 * (i + j)].push_back(x0);
            // add x start (x0) to trajectories
            s.trajectories[2 * (i + j) + 1].push_back(y0);
            // add y start (y0) to trajectories
        }
        double t = 0;  // t start
        if constexpr (Hybrid) {
            // This feature should only really be used with more accurate
            // numerical methods like rungekutta4
            double dtFast = 10 * dt_min;
            // Use larger time step
            while (stopFunctionHybrid(variables, dtFast)) {
                /*numericalMethod<Numerical, singleTrajDims>(variables,
                                                            dtFast, k, cw_2);
                t += dtFast; // adjust time
                if constexpr (AddTraj) {
                    s.trajectories[2 * (i +
                j)].push_back(variables[singleTrajVars::x]); s.trajectories[2 *
                (i + j) + 1].push_back(variables[singleTrajVars::y]);
                }*/
                singleTrajStep<AddTraj, Numerical, singleTrajDims>(
                    &shellCalc::stopFunctionHybrid, s.trajectories[2 * (i + j)],
                    s.trajectories[2 * (i + j) + 1], variables, t, dt_min, k,
                    cw_2);
            }
        }

        while (stopFunction(variables, dt_min)) {
            singleTrajStep<AddTraj, Numerical, singleTrajDims>(
                &shellCalc::stopFunction, s.trajectories[2 * (i + j)],
                s.trajectories[2 * (i + j) + 1], variables, t, dt_min, k, cw_2);
        }
        s.get_impact(i + j, impact::distance) = variables[singleTrajVars::x];
        vx[j] = variables[singleTrajVars::v_x];
        vy[j] = variables[singleTrajVars::v_y];
        tVec[j] = t;
    }

    template <bool AddTraj, unsigned int Numerical>
    void multiTraj(const unsigned int i, const unsigned int j, shell &s,
                   double *vx, double *vy, double *tVec) {
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();

        std::array<double, vSize> groupX;
        std::array<double, vSize> groupY;
        std::array<double, vSize> groupVX;
        std::array<double, vSize> groupVY;

        auto checkContinue = [&]() {
            bool any = false;
            for (unsigned int i = 0; i < vSize; i++) {
                bool test = groupY[i] >= 0;
                if (test) {
                    any = test;
                }
            }
            return any;
        }

        while(checkContinue){
            for(unsigned int i=0; i< vSize; i++){
                numericalMethod<Numerical, 2>({groupX[i], groupY[i], groupVX[i], groupVY[i]}, dt_min, k, cw_2);
            }
        }
    }

    // Several trajectories done in one chunk to allow for vectorization
    template <bool AddTraj, unsigned int Numerical, bool Hybrid, bool Fit,
              bool nonAP>
    void impactGroup(const unsigned int i, shell *const shellPointer) {
        shell &s = *shellPointer;
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_normalizationR();

        double vx[vSize], vy[vSize], tVec[vSize];
        for (unsigned int j = 0; j < vSize; j++) {
            if constexpr (!Fit) {
                s.get_impact(i + j, impact::launchAngle) =
                    std::fma(precision, (i + j), min);
            }
            double radianLaunch =
                s.get_impact(i + j, impact::launchAngle) * M_PI / 180;
            vx[j] = s.get_v0() * cos(radianLaunch);
            vy[j] = s.get_v0() * sin(radianLaunch);
        }
        for (unsigned int j = 0; (j + i < s.impactSize) & (j < vSize); j++) {
            singleTraj<AddTraj, Numerical, Hybrid>(i, j, s, vx, vy, tVec);
        }
        for (unsigned int j = 0; j < vSize; j++) {
            double IA_R = atan(vy[j] / vx[j]);

            s.get_impact(i + j, impact::impactAngleHorizontalRadians) = IA_R;
            double IAD_R = M_PI_2 + IA_R;
            double IA_D = IA_R * 180 / M_PI;
            s.get_impact(i + j, impact::impactAngleHorizontalDegrees) =
                IA_D * -1;
            s.get_impact(i + j, impact::impactAngleDeckDegrees) = 90 + IA_D;

            double IV = sqrt(vx[j] * vx[j] + vy[j] * vy[j]);
            s.get_impact(i + j, impact::impactVelocity) = IV;
            s.get_impact(i + j, impact::timeToTarget) = tVec[j];
            s.get_impact(i + j, impact::timeToTargetAdjusted) = tVec[j] / 3.1;

            if constexpr (!Fit) {
                if constexpr (nonAP) {
                    s.get_impact(i + j, impact::rawPenetration) = s.nonAP;
                    s.get_impact(i + j,
                                 impact::effectivePenetrationHorizontal) =
                        s.nonAP;
                    s.get_impact(i + j, impact::effectivePenetrationDeck) =
                        s.nonAP;
                    s.get_impact(
                        i + j,
                        impact::effectivePenetrationHorizontalNormalized) =
                        s.nonAP;
                    s.get_impact(i + j,
                                 impact::effectivePenetrationDeckNormalized) =
                        s.nonAP;
                } else {
                    double rawPenetration = pPPC * pow(IV, 1.1001);
                    s.get_impact(i + j, impact::rawPenetration) =
                        rawPenetration;

                    s.get_impact(i + j,
                                 impact::effectivePenetrationHorizontal) =
                        rawPenetration * cos(IA_R);
                    s.get_impact(i + j, impact::effectivePenetrationDeck) =
                        rawPenetration * cos(IAD_R);

                    s.get_impact(
                        i + j,
                        impact::effectivePenetrationHorizontalNormalized) =
                        rawPenetration *
                        cos(calcNormalizationR(IA_R, normalizationR));
                    s.get_impact(i + j,
                                 impact::effectivePenetrationDeckNormalized) =
                        rawPenetration *
                        cos(calcNormalizationR(IAD_R, normalizationR));
                }
            }
        }
    }

   public:
    unsigned int calculateAlignmentSize(unsigned int unalignedSize) {
        return vSize - (unalignedSize % vSize) + unalignedSize;
    }
    // Templates to reduce branching
    template <unsigned int Numerical, bool Hybrid>
    void calculateImpact(
        shell &s, bool addTraj,
        unsigned int nThreads = std::thread::hardware_concurrency()) {
        if (addTraj) {
            calculateImpact<true, Numerical, Hybrid>(s, nThreads);
        } else {
            calculateImpact<false, Numerical, Hybrid>(s, nThreads);
        }
    }

    template <bool AddTraj, unsigned int Numerical, bool Hybrid>
    void calculateImpact(
        shell &s, unsigned int nThreads = std::thread::hardware_concurrency()) {
        if (s.enableNonAP) {
            calculateImpact<true, Numerical, Hybrid, true>(s, nThreads);
        } else {
            calculateImpact<true, Numerical, Hybrid, false>(s, nThreads);
        }
    }

    template <bool AddTraj, unsigned int Numerical, bool Hybrid, bool nonAP>
    void calculateImpact(
        shell &s, unsigned int nThreads = std::thread::hardware_concurrency()) {
        s.impactSize = ((max / precision - min / precision)) + 1;
        s.impactSizeAligned = calculateAlignmentSize(s.impactSize);
        if constexpr (AddTraj) {
            s.trajectories.resize(2 * s.impactSize);
        }
        s.impactData.resize(impact::maxColumns * s.impactSizeAligned);

        if (nThreads > std::thread::hardware_concurrency()) {
            nThreads = std::thread::hardware_concurrency();
        }
        unsigned int length = ceil((double)s.impactSize / vSize);
        unsigned int assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(
            assigned, length, s.impactSize, this,
            &shellCalc::impactGroup<AddTraj, Numerical, Hybrid, false, nonAP>,
            &s);

        s.completedImpact = true;
    }

    template <unsigned int Numerical>
    void calculateFit(
        shell &s, unsigned int nThreads = std::thread::hardware_concurrency()) {
        if (nThreads > std::thread::hardware_concurrency()) {
            nThreads = std::thread::hardware_concurrency();
        }
        unsigned int length = ceil((double)s.impactSize / vSize);
        unsigned int assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(assigned, length, s.impactSize, this,
                         &shellCalc::multiTraj<false, Numerical, false, true>,
                         &s);
    }

   private:
    void checkRunImpact(shell &s) {
        if (!s.completedImpact) {
            std::cout << "Standard Not Calculated - Running automatically\n";
            calculateImpact<false, numerical::forwardEuler, false>(s);
        }
    }

    // Check Angles Section
    // template <short fusing> explanation:
    // Fusing is done using templates to reduce in loop branching and
    // computational time in some cases.

    // Possible Values: 0 - Never Fusing 1 - Check 2 - Always Fusing
    enum fuseStatus { never, check, always };
    template <short fusing, bool nonAP, bool nonAPPerforated,
              bool disableRicochet>
    void multiAngles(const unsigned int i, const double thickness,
                     const double inclination_R, const double fusingAngle,
                     shell *const shellPointer) {
        static_assert(fusing <= 2 && fusing >= 0, "Invalid fusing parameter");
        shell &s = *shellPointer;
        for (unsigned int j = 0; j < vSize; j++) {
            double fallAngleAdjusted =
                s.get_impact(
                    i + j,
                    impact::impactDataIndex::impactAngleHorizontalRadians) +
                inclination_R;
            double rawPenetration =
                s.get_impact(i + j, impact::impactDataIndex::rawPenetration);

            double penetrationCriticalAngle;
            if constexpr (nonAP) {
                if constexpr (nonAPPerforated) {
                    penetrationCriticalAngle = M_PI_2;
                } else {
                    penetrationCriticalAngle = 0;
                }
            } else {
                penetrationCriticalAngle =

                    (acos(thickness / rawPenetration) + s.get_normalizationR());

                penetrationCriticalAngle = std::isnan(penetrationCriticalAngle)
                                               ? 0
                                               : penetrationCriticalAngle;
            }

            std::array<double, 4> criticalAngles;
            if constexpr (disableRicochet) {
                criticalAngles = {M_PI_2, M_PI_2, penetrationCriticalAngle,
                                  fusingAngle};
            } else {
                criticalAngles = {s.ricochet0R, s.ricochet1R,

                                  penetrationCriticalAngle, fusingAngle};
            }
            std::array<double, 4> out;
            for (unsigned int k = 0; k < 2; k++) {
                out[k] = acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                out[k] = std::isnan(out[k]) ? 0 : out[k];
            }

            {
                int k = angle::armorRadians / 2;

                if (criticalAngles[k] < M_PI_2) {
                    out[k] =
                        acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                    out[k] = std::isnan(out[k]) ? 0 : out[k];
                } else {
                    out[k] = M_PI_2;
                }
            }
            {
                int k = angle::fuseRadians / 2;
                if constexpr (fusing == fuseStatus::never) {
                    out[k] = M_PI_2;
                } else if (fusing == fuseStatus::check) {
                    out[k] =
                        acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                    out[k] = std::isnan(out[k]) ? 0 : out[k];
                } else if (fusing == fuseStatus::always) {
                    out[k] = 0;
                }
            }

            s.get_angle(i + j, angle::angleDataIndex::ricochetAngle0Radians) =
                out[0];
            s.get_angle(i + j, angle::angleDataIndex::ricochetAngle1Radians) =
                out[1];
            s.get_angle(i + j, angle::angleDataIndex::armorRadians) = out[2];
            s.get_angle(i + j, angle::angleDataIndex::fuseRadians) = out[3];

            for (unsigned int k = 0; k < angle::maxColumns / 2; k++) {
                out[k] *= 180 / M_PI;
            }

            s.get_angle(i + j, angle::angleDataIndex::ricochetAngle0Degrees) =
                out[0];
            s.get_angle(i + j, angle::angleDataIndex::ricochetAngle1Degrees) =
                out[1];
            s.get_angle(i + j, angle::angleDataIndex::armorDegrees) = out[2];
            s.get_angle(i + j, angle::angleDataIndex::fuseDegrees) = out[3];
        }
    }

   public:
    void calculateAngles(
        const double thickness, const double inclination, shell &s,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        if (s.enableNonAP) {
            if (s.nonAP >= thickness) {
                calculateAngles<true, true>(thickness, inclination, s,
                                            nThreads);
            } else {
                calculateAngles<true, false>(thickness, inclination, s,
                                             nThreads);
            }
        } else {
            calculateAngles<false, false>(thickness, inclination, s, nThreads);
        }
    }
    template <bool nonAP, bool nonAPPerforated>
    void calculateAngles(
        const double thickness, const double inclination, shell &s,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        if (s.ricochet0 >= 90) {
            calculateAngles<nonAP, nonAPPerforated, true>(
                thickness, inclination, s, nThreads);
        } else {
            calculateAngles<nonAP, nonAPPerforated, false>(
                thickness, inclination, s, nThreads);
        }
    }

    template <bool nonAP, bool nonAPPerforated, bool disableRicochet>
    void calculateAngles(
        const double thickness, const double inclination, shell &s,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        checkRunImpact(s);

        s.angleData.resize(angle::maxColumns * s.impactSizeAligned);
        std::copy_n(s.get_impactPtr(0, impact::distance), s.impactSize,
                    s.get_anglePtr(0, angle::distance));

        unsigned int length = (unsigned int)ceil((double)s.impactSize / vSize);
        unsigned int assigned = assignThreadNum(length, nThreads);

        double inclination_R = inclination / 180 * M_PI;
        double fusingAngle;
        fusingAngle = acos(thickness / s.threshold) + s.get_normalizationR();
        if (std::isnan(fusingAngle)) {
            mtFunctionRunner(
                assigned, length, s.impactSize, this,
                &shellCalc::multiAngles<fuseStatus::always, nonAP,
                                        nonAPPerforated, disableRicochet>,
                thickness, inclination_R, fusingAngle, &s);
        } else {
            if (fusingAngle > M_PI_2) {
                mtFunctionRunner(
                    assigned, length, s.impactSize, this,
                    &shellCalc::multiAngles<fuseStatus::never, nonAP,
                                            nonAPPerforated, disableRicochet>,
                    thickness, inclination_R, fusingAngle, &s);
            } else {
                mtFunctionRunner(
                    assigned, length, s.impactSize, this,
                    &shellCalc::multiAngles<fuseStatus::check, nonAP,
                                            nonAPPerforated, disableRicochet>,
                    thickness, inclination_R, fusingAngle, &s);
            }
        }
        s.completedAngles = true;
    }

    // Post-Penetration Section

   private:
    template <bool fast>
    void postPenTraj(const unsigned int i, shell &s, double v_x, double v_y,
                     double v_z, double thickness) {
        const double notFusedCode = -1;
        if constexpr (fast) {
            bool positive = v_x > 0;
            double x = v_x * s.fuseTime * positive;
            s.get_postPen(i, post::x, 0) = x;
            s.get_postPen(i, post::y, 0) = v_y * s.fuseTime * positive;
            s.get_postPen(i, post::z, 0) = v_z * s.fuseTime * positive;

            bool fuse = thickness >= s.threshold;
            s.get_postPen(i, post::xwf, 0) = (fuse)*x + !(fuse)*notFusedCode;
        } else {
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
            double pos[3], velocities[3], velocitiesSquared[3],
                dragIntermediary[3];
            double xz_dragIntermediary[2];
            pos[0] = xf0, pos[1] = yf0, pos[2] = xf0;
            velocities[0] = v_x, velocities[1] = v_y, velocities[2] = v_z;
            t = 0;
            if (v_x > 0) {
                while (t < s.fuseTime) {
                    for (int l = 0; l < 3; l++) {
                        pos[l] += velocities[l] * dtf;
                    }
                    // Calculate air density - likely unnecessary for this
                    // section as distances are so short
                    T = t0 - L * pos[1];
                    p = p0 * pow((1 - L * pos[1] / t0), gMRL);
                    rho = p * M / (R * T);

                    // Calculated drag deceleration

                    for (int l = 0; l < 3; l++) {
                        velocitiesSquared[l] = velocities[l] * velocities[l];
                    }
                    // velocitiesSquared = velocities * velocities
                    xz_dragIntermediary[0] =
                        (k * rho) *
                        (cw_1 * velocitiesSquared[0] + cw_2 * velocities[0]);
                    xz_dragIntermediary[1] =
                        (k * rho) *
                        (cw_1 * velocitiesSquared[2] + cw_2 * velocities[2]);
                    // xz_dragIntermediary = (k * rho) * (cw_1 *
                    // velocitiesSquared[2, 0] + cw_2 * velocities[2, 0])
                    dragIntermediary[0] = xz_dragIntermediary[0];  // x
                    dragIntermediary[1] =
                        (g + k * rho *
                                 (cw_1 * velocitiesSquared[1] +
                                  cw_2 * fabs(velocities[1])) *
                                 signum(velocities[1]));
                    dragIntermediary[2] = xz_dragIntermediary[1];  // z
                    for (int l = 0; l < 3; l++) {
                        velocities[l] -= dtf * dragIntermediary[l];
                    }
                    // velocities -= dtf * dragIntermediary
                    t += dtf;
                }
                s.get_postPen(i, post::x, 0) = pos[0];
                s.get_postPen(i, post::y, 0) = pos[1];
                s.get_postPen(i, post::z, 0) = pos[2];
                s.get_postPen(i, post::xwf, 0) =
                    (thickness >= s.threshold) * pos[0] +
                    !(thickness >= s.threshold) * notFusedCode;
            } else {
                s.get_postPen(i, post::x, 0) = 0;
                s.get_postPen(i, post::y, 0) = 0;
                s.get_postPen(i, post::z, 0) = 0;
                s.get_postPen(i, post::xwf, 0) = 0;
            }
        }
    }

    template <bool changeDirection, bool fast>
    void multiPostPen(unsigned int i, const double thickness,
                      const double inclination_R, shell *const shellPointer) {
        shell &s = *shellPointer;
        double hAngleV[vSize], vAngleV[vSize];
        double v0V[vSize], penetrationV[vSize], eThicknessV[vSize];
        double v_x[vSize], v_y[vSize], v_z[vSize];
        unsigned int distIndex = (i < s.impactSize) ? i : i % s.impactSize;
        unsigned int anglesIndex = i / s.impactSize;

        unsigned int j, k = 0;

        if (i + vSize <= s.postPenSize) {
            std::copy_n(s.get_postPenPtr(i, post::angle, 0), vSize, hAngleV);
        } else {
            for (j = 0; (i + j) < s.postPenSize; j++) {
                hAngleV[j] = s.postPenData[i + j];
            }
        }

        if (distIndex < s.impactSize - vSize + 1) {
            std::copy_n(s.get_impactPtr(distIndex,
                                        impact::impactAngleHorizontalRadians),
                        vSize, vAngleV);
            std::copy_n(s.get_impactPtr(distIndex, impact::rawPenetration),
                        vSize, penetrationV);
            std::copy_n(s.get_impactPtr(distIndex, impact::impactVelocity),
                        vSize, v0V);
        } else {
            for (j = 0; (j + distIndex < s.impactSize) && (j < vSize); j++) {
                vAngleV[j] = s.get_impact(distIndex + j,
                                          impact::impactAngleHorizontalRadians);
                penetrationV[j] =
                    s.get_impact(distIndex + j, impact::rawPenetration);
                v0V[j] = s.get_impact(distIndex + j, impact::impactVelocity);
            }
            if (anglesIndex < s.postPenSize / s.impactSize) {
                for (; (j < vSize); j++) {
                    vAngleV[j] =
                        s.get_impact(k, impact::impactAngleHorizontalRadians);
                    penetrationV[j] = s.get_impact(k, impact::rawPenetration);
                    v0V[j] = s.get_impact(k, impact::impactVelocity);
                    k++;
                }
            }
        }

        for (unsigned int l = 0; l < vSize; l++) {
            double HA_R = hAngleV[l] * M_PI / 180;     // lateral  angle radians
            double VA_R = vAngleV[l] + inclination_R;  // vertical angle radians
            double cAngle = acos(cos(HA_R) * cos(VA_R));
            double nCAngle = calcNormalizationR(cAngle, s.get_normalizationR());
            double eThickness = thickness / cos(nCAngle);
            double pPV = v0V[l] * (1 - exp(1 - penetrationV[l] / eThickness));

            if constexpr (changeDirection) {
                double hFAngle = atan(tan(nCAngle) * tan(HA_R) / tan(cAngle));
                double vFAngle = atan(tan(nCAngle) * cos(hFAngle) * tan(VA_R) /
                                      cos(HA_R) / tan(cAngle));

                double v_x0 = pPV * cos(vFAngle) * cos(hFAngle);
                double v_y0 = pPV * cos(vFAngle) * sin(hFAngle);

                v_x[l] = v_x0 * cos(inclination_R) + v_y0 * sin(inclination_R);
                v_z[l] = v_y0 * cos(inclination_R) + v_x0 * sin(inclination_R);
                v_y[l] = pPV * sin(vFAngle);
            } else {
                v_x[l] = pPV * cos(VA_R) * cos(HA_R);
                v_z[l] = pPV * cos(VA_R) * sin(HA_R);
                v_y[l] = pPV * sin(VA_R);
            }
            eThicknessV[l] = eThickness;
        }

        for (unsigned int j = 0; (j < vSize) & (j + i < s.postPenSize); j++) {
            postPenTraj<fast>(i + j, s, v_x[j], v_y[j], v_z[j], eThicknessV[j]);
        }
        // std::cout<<index<<" Completed\n";
    }

    // Probably unnecessary...
    void fillCopy(std::atomic<int> &counter, unsigned int assigned,
                  unsigned int id, shell *const s,
                  std::vector<double> *angles) {
        for (unsigned int i = angles->size() * id / assigned;
             i < angles->size() * (id + 1) / assigned; i++) {
            std::fill_n(s->get_postPenPtr(0, post::angle, i), s->impactSize,
                        (double)(*angles)[i]);
            std::copy_n(
                s->get_impactPtr(0, impact::distance), s->impactSize,
                s->postPenData.begin() + s->postPenSize + i * s->impactSize);
        }
        counter.fetch_add(1, std::memory_order_relaxed);
    }

    void parallelFillCopy(shell *const s, std::vector<double> *angles,
                          unsigned int nThreads) {
        std::vector<std::thread> threads;
        std::atomic<int> counter{0};
        unsigned int length = angles->size();
        unsigned int assigned = assignThreadNum(length, nThreads);
        for (unsigned int i = 0; i < assigned - 1; i++) {
            threads.push_back(std::thread(
                [=, &counter] { fillCopy(counter, assigned, i, s, angles); }));
        }
        fillCopy(counter, assigned, assigned - 1, s, angles);
        for (unsigned int i = 0; i < assigned - 1; i++) {
            threads[i].join();
        }
    }

   public:
    // Again templates to reduce branching
    void calculatePostPen(
        const double thickness, const double inclination, shell &s,
        std::vector<double> &angles, const bool changeDirection = true,
        const bool fast = false,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        // Specifies whether normalization alters the trajectory of the
        // shell Though effect is not too significant either way
        if (changeDirection) {
            calculatePostPen<true>(thickness, inclination, s, angles, fast,
                                   nThreads);
        } else {
            calculatePostPen<false>(thickness, inclination, s, angles, fast,
                                    nThreads);
        }
    }

    template <bool changeDirection>
    void calculatePostPen(
        const double thickness, const double inclination, shell &s,
        std::vector<double> &angles, const bool fast = false,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        /* Specifies whether to perform ballistic calculations or just
           multiply velocity by fusetime.
           It is faster but post-penetration already runs really fast... */
        if (fast) {
            calculatePostPen<changeDirection, true>(thickness, inclination, s,
                                                    angles, nThreads);
        } else {
            calculatePostPen<changeDirection, false>(thickness, inclination, s,
                                                     angles, nThreads);
        }
    }

    template <bool changeDirection, bool fast>
    void calculatePostPen(const double thickness, const double inclination,
                          shell &s, std::vector<double> &angles,
                          const unsigned int nThreads) {
        checkRunImpact(s);

        s.postPenSize = s.impactSize * angles.size();
        s.postPenData.resize(6 * s.postPenSize);

        parallelFillCopy(&s, &angles, nThreads);
        double inclination_R = M_PI / 180 * inclination;
        unsigned int length = ceil((double)s.postPenSize / vSize);
        unsigned int assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(assigned, length, s.postPenSize, this,
                         &shellCalc::multiPostPen<changeDirection, fast>,
                         thickness, inclination_R, &s);

        s.completedPostPen = true;
    }
};

}  // namespace shell