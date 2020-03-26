#define _USE_MATH_DEFINES
#include <cmath>

#include "concurrentqueue/concurrentqueue.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace shell {

namespace impact {
static constexpr unsigned int maxColumns = 13;
enum impactDataIndex {
    distance,
    launchA,
    impactAHR,
    impactAHD,
    impactV,
    rawPen,
    ePenH,
    ePenHN,
    impactADD,
    ePenD,
    ePenDN,
    tToTarget,
    tToTargetA
};
static_assert(tToTargetA == (maxColumns - 1), "Invalid standard columns");
} // namespace impact

namespace angle {
static constexpr unsigned int maxColumns = 9;
enum angleDataIndex {
    distance,
    ra0,
    ra0D,
    ra1,
    ra1D,
    armor,
    armorD,
    fuse,
    fuseD
};
static_assert(fuseD == (maxColumns - 1), "Invalid angle columns");
} // namespace angle

namespace post {
static constexpr unsigned int maxColumns = 6;
enum postPenDataIndex { angle, distance, x, y, z, xwf };
static_assert(xwf == (maxColumns - 1), "Invaild postpen columns");
} // namespace post

enum numerical { forwardEuler, rungeKutta4 };

class shell {
private:                  // Description         units
    double v0;            // muzzle velocity     m/s
    double caliber;       // shell caliber       m
    double krupp;         // shell krupp         [ndim]
    double mass;          // shell mass          kg
    double cD;            // coefficient of drag [ndim]
    double normalization; // shell normalization degrees
    double threshold;     // fusing threshold    mm
    double fuseTime;      // fuse time           s
    double ricochet0;     // ricochet angle 0    degrees
    double ricochet1;     // ricochet angle 1    degrees
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
        pPPC = 0.5561613 * krupp / 2400 * pow(mass, 0.55) /
               pow((caliber * 1000), 0.65);
        // condensed penetration coefficient
        normalizationR = normalization / 180 * M_PI;
        // normalization (radians)
        ricochet0R = ricochet0 / 180 * M_PI;
        ricochet1R = ricochet1 / 180 * M_PI;
    }

    unsigned int impactSize,
        postPenSize; // number of distances in: standard, postPen
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
          const double ricochet1, const std::string &name) {
        setValues(caliber, v0, cD, mass, krupp, normalization, fuseTime,
                  threshold, ricochet0, ricochet1, name);
    }

    void setValues(const double caliber, const double v0, const double cD,
                   const double mass, const double krupp,
                   const double normalization, const double fuseTime,
                   const double threshold, const double ricochet0,
                   const double ricochet1, const std::string &name) {
        //                                                    Impact   PostPen
        this->fuseTime = fuseTime; // Shell fusetime        | No     | Yes
        this->v0 = v0;             // Shell muzzle velocity | Yes    | No
        this->caliber = caliber;   // Shell caliber         | Yes    | Yes
        this->krupp = krupp;       // Shell krupp           | Yes    | Yes
        this->mass = mass;         // Shell mass            | Yes    | Yes
        this->normalization =
            normalization; // Shell normalization   | Yes    | Yes
        this->cD = cD;     // Shell air drag coefficient | Yes    | Yes
        this->name = name; // Shell name  | No     | No
        this->ricochet0 = ricochet0; // Ricochet Angle 0 | No     | Yes
        this->ricochet1 = ricochet1; // Ricochet Angle 1 | No     | Yes
                                     // Shell fusing threshold | No     | Yes
        this->threshold = threshold;

        preProcess();
    }
    // Setter Functions
    // Note: Be sure to call preprocess after changing to make sure calculations
    // are correct
    void set_v0(const double v0) { this->v0 = v0; }
    void set_caliber(const double caliber) { this->caliber = caliber; }
    void set_krupp(const double krupp) { this->krupp = krupp; }
    void set_mass(const double mass) { this->mass = mass; }
    void set_normalization(const double normalization) {
        this->normalization = normalization;
    }
    void set_cD(const double cD) { this->cD = cD; }
    void set_name(const std::string &name) { this->name = name; }
    void set_threshold(const double threshold) { this->threshold = threshold; }
    void set_fuseTime(const double fuseTime) { this->fuseTime = fuseTime; }
    void set_ricochet0(const double ricochet0) { this->ricochet0 = ricochet0; }
    void set_ricochet1(const double ricochet1) { this->ricochet1 = ricochet1; }

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

    // fixed
    const double get_v0() { return v0; }
    const double get_k() { return k; }
    const double get_cw_2() { return cw_2; }
    const double get_pPPC() { return pPPC; }
    const double get_normalizationR() { return normalizationR; }

    // changeable
    double &get_threshold() { return threshold; }
    double &get_fuseTime() { return fuseTime; }
    double &get_ricochet0() { return ricochet0; }
    double &get_ricochet1() { return ricochet1; }
    double &get_ricochet0R() { return ricochet0R; }
    double &get_ricochet1R() { return ricochet1R; }

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
    double g = 9.81;      // Gravitational Constant       | m/(s^2)
    double t0 = 288;      // Temperature at Sea Level     | K
    double L = 0.0065;    // Atmospheric Lapse Rate       | C/m
    double p0 = 101325;   // Pressure at Sea Level        | Pa
    double R = 8.31447;   // Ideal Gas Constant           | J/(mol K)
    double M = 0.0289644; // Molarity of Air at Sea Level | kg/mol
    double cw_1 = 1;
    // Calculation Parameters
    double max = 25;       // Max Angle                    | degrees
    double min = 0;        // Min Angle                    | degrees
    double precision = .1; // Angle Step                   | degrees
    double x0 = 0, y0 = 0; // Starting x0, y0              | m
    double dt_min = .01;   // Time step                    | s

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
        "Size of double is not 8 - required for AVX2"); // Use float64
                                                        // in the future
    static_assert(std::numeric_limits<double>::is_iec559,
                  "Type is not IEE754 compliant");
    static constexpr unsigned int vSize = (256 / 8) / sizeof(double);

public:
    double calcNormalizationR(const double angle,
                              const double normalizationR) { // Input in radians
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
    // mini 'threadpool' used to kick off multithreaded functions
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
            moodycamel::ConcurrentQueue<int> workQueue;
            std::vector<std::thread> threads(assigned - 1);
            for (unsigned int i = 0; i < assigned - 1; i++) {
                threads[i] = std::thread([&] {
                    mtWorker(counter, workQueue, length, i, object, function,
                             args...);
                });
            }

            int buffer[workQueueBufferSize];
            unsigned int bCounter = 0;
            for (unsigned int i = 0; i < size; i += vSize) {
                buffer[bCounter] = i;
                bCounter++;
                if (bCounter == workQueueBufferSize) {
                    workQueue.enqueue_bulk(buffer, bCounter);
                    bCounter = 0;
                }
            }
            workQueue.enqueue_bulk(buffer, bCounter);

            mtWorker(counter, workQueue, length, assigned - 1, object, function,
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
    void mtWorker(std::atomic<int> &counter,
                  moodycamel::ConcurrentQueue<int> &workQueue, const int length,
                  const int threadID, O object, F function, Args... args) {
        // threadID is largely there for debugging
        while (counter < length) {
            int index;
            if (workQueue.try_dequeue(index)) {
                (object->*function)(index, args...);
                counter.fetch_add(1, std::memory_order_relaxed);
            } else {
                std::this_thread::yield();
            }
        }
    }

    // 0: x 1: y 2: v_x 3: v_y
    template <unsigned int dims>
    std::array<double, dims * 2>
    calcDeltas(std::array<double, dims * 2> current, const double dt,
               const double k, const double cw_2) {
        enum position { x, y, z };
        enum velocity { v_x, v_y, v_z };

        std::array<double, dims * 2> deltas;
        for (unsigned int i = 0; i < dims; i++) {
            deltas[i] = current[i + dims] * dt;
            current[i] += deltas[i];
        }

        double T, p, rho;
        // Calculating air density
        T = t0 - L * current[position::y];
        // Calculating air temperature at altitude
        p = p0 * pow((1 - L * current[position::y] / t0), (g * M / (R * L)));
        // Calculating air pressure at altitude
        rho = p * M / (R * T);
        // Use ideal gas law to calculate air density

        // Calculate drag deceleration
        std::array<double, dims> velocitySquared;
        for (unsigned int l = 0; l < dims; l++) {
            velocitySquared[l] = current[dims + l] * current[dims + l];
        } // v^2 = v * v

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

        deltas[dims + velocity::v_y] = g - deltas[dims + velocity::v_y];
        for (unsigned int i = 0; i < dims; i++) { // v -= drag * dt
            deltas[dims + i] *= (-1 * dt);
        }
        return deltas;
    }
    template <std::size_t N>
    void fmaArr(double x, std::array<double, N> y, std::array<double, N> &z) {
        for (unsigned int i = 0; i < N; i++) {
            z[i] = std::fma(x, y[i], z[i]);
        }
    }

    // Numerical Methods
    template <unsigned int dims>
    void forwardEuler(std::array<double, 2 * dims> &input, const double dt,
                      const double k, const double cw_2) {
        std::array<double, 2 *dims> delta =
            calcDeltas<singleTrajDims>(input, dt, k, cw_2);
        for (unsigned int i = 0; i < 2 * dims; i++) {
            input[i] += delta[i];
        }
    }

    template <unsigned int dims>
    void rungeKutta4(std::array<double, 2 * dims> &input, const double dt,
                     const double k, const double cw_2) {
        std::array<double, 2 * dims> k1, k2, k3, k4;
        std::array<double, 2 *dims> intermediate = input;
        k1 = calcDeltas<dims>(input, dt, k, cw_2);
        fmaArr<2 * dims>(.5, k1, intermediate);
        k2 = calcDeltas<dims>(intermediate, dt, k, cw_2);
        intermediate = input;
        fmaArr<2 * dims>(.5, k2, intermediate);
        k3 = calcDeltas<dims>(intermediate, dt, k, cw_2);
        intermediate = input;
        fmaArr<2 * dims>(1, k3, intermediate);
        k4 = calcDeltas<dims>(intermediate, dt, k, cw_2);

        fmaArr<2 * dims>(2, k2, k1);
        fmaArr<2 * dims>(2, k3, k1);
        fmaArr<2 * dims>(1, k4, k1);
        fmaArr<2 * dims>((1.0 / 6.0), k1, input);
    }

    // Please edit numerical enum before changing this
    template <unsigned int Numerical, unsigned int Dims>
    void numericalMethod(std::array<double, 4> &input, const double dt,
                         const double k, const double cw_2) {
        if constexpr (Numerical == numerical::forwardEuler) {
            forwardEuler<Dims>(input, dt, k, cw_2);
        } else if (Numerical == numerical::rungeKutta4) {
            rungeKutta4<Dims>(input, dt, k, cw_2);
        }
    }
    // Impact Calculations Region
    enum singleTrajVars { x, y, v_x, v_y };
    static constexpr unsigned int singleTrajDims = 2;

    template <bool AddTraj, unsigned int Numerical, bool Hybrid>
    void singleTraj(const unsigned int i, const unsigned int j, shell &s,
                    double *vx, double *vy, double *tVec) {
        static constexpr unsigned int __TrajBuffer__ = 128;
        // likely unnecessary - just part of a feature that was never
        // implemented
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();

        unsigned int counter = 0;
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
        std::array<double, __TrajBuffer__> xT, yT;
        // setting initial values
        std::array<double, 4> variables{x0, y0, vx[j], vy[j]};
        s.trajectories[2 * (i + j)].push_back(x0);
        // add x start (x0) to trajectories
        s.trajectories[2 * (i + j) + 1].push_back(y0);
        // add y start (y0) to trajectories
        double t = 0; // t start
        if constexpr (Hybrid) {
            // This feature should only really be used with more accurate
            // numerical methods like rungekutta4
            double dtFast = 8 * dt_min;
            // Use larger time step
            while (variables[singleTrajVars::y] +
                       variables[singleTrajVars::v_y] * dtFast -
                       g / 2 * dtFast * dtFast >=
                   0) {
                for (counter = 0; (counter < __TrajBuffer__) &
                                  (variables[singleTrajVars::y] >= 0);
                     counter++) {

                    numericalMethod<Numerical, singleTrajDims>(variables,
                                                               dtFast, k, cw_2);
                    t += dtFast; // adjust time

                    if constexpr (AddTraj) {
                        xT[counter] = variables[singleTrajVars::x];
                        yT[counter] = variables[singleTrajVars::y];
                    }
                }
                if constexpr (AddTraj) {
                    s.trajectories[2 * (i + j)].insert(
                        s.trajectories[2 * (i + j)].end(), xT.data(),
                        &xT[counter]);
                    s.trajectories[2 * (i + j) + 1].insert(
                        s.trajectories[2 * (i + j) + 1].end(), yT.data(),
                        &yT[counter]);
                }
            }
        }

        while (variables[singleTrajVars::y] >= 0) {
            for (counter = 0; (counter < __TrajBuffer__) &
                              (variables[singleTrajVars::y] >= 0);
                 counter++) {

                numericalMethod<Numerical, singleTrajDims>(variables, dt_min, k,
                                                           cw_2);
                t += dt_min; // adjust time

                if constexpr (AddTraj) {
                    xT[counter] = variables[singleTrajVars::x];
                    yT[counter] = variables[singleTrajVars::y];
                }
            }
            if constexpr (AddTraj) {
                s.trajectories[2 * (i + j)].insert(
                    s.trajectories[2 * (i + j)].end(), xT.data(), &xT[counter]);
                s.trajectories[2 * (i + j) + 1].insert(
                    s.trajectories[2 * (i + j) + 1].end(), yT.data(),
                    &yT[counter]);
            }
        }
        s.get_impact(i + j, impact::distance) = variables[singleTrajVars::x];
        vx[j] = variables[singleTrajVars::v_x];
        vy[j] = variables[singleTrajVars::v_y];
        tVec[j] = t;
    }

    // Several trajectories done in one chunk to allow for vectorization
    template <bool AddTraj, unsigned int Numerical, bool Hybrid>
    void multiTraj(const unsigned int i, shell *const shellPointer) {
        shell &s = *shellPointer;
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_normalizationR();

        double vx[vSize], vy[vSize], tVec[vSize];
        for (unsigned int j = 0; j < vSize; j++) {
            s.get_impact(i + j, impact::launchA) =
                std::fma(precision, (i + j), min);
            double radianLaunch =
                s.get_impact(i + j, impact::launchA) * M_PI / 180;
            vx[j] = s.get_v0() * cos(radianLaunch);
            vy[j] = s.get_v0() * sin(radianLaunch);
        }

        for (unsigned int j = 0; (j + i < s.impactSize) & (j < vSize); j++) {
            singleTraj<AddTraj, Numerical, Hybrid>(i, j, s, vx, vy, tVec);
        }

        for (unsigned int j = 0; j < vSize; j++) {
            double IA_R = atan(vy[j] / vx[j]);
            s.get_impact(i + j, impact::impactAHR) = IA_R;
            double IAD_R = M_PI / 2 + IA_R;
            double IA_D = IA_R * 180 / M_PI;
            s.get_impact(i + j, impact::impactAHD) = IA_D;
            s.get_impact(i + j, impact::impactADD) = 90 + IA_D;

            double IV = sqrt(vx[j] * vx[j] + vy[j] * vy[j]);
            s.get_impact(i + j, impact::impactV) = IV;
            double rawPen = pPPC * pow(IV, 1.1);
            s.get_impact(i + j, impact::rawPen) = rawPen;

            s.get_impact(i + j, impact::ePenH) = rawPen * cos(IA_R);
            s.get_impact(i + j, impact::ePenD) = rawPen * cos(IAD_R);

            s.get_impact(i + j, impact::ePenHN) =
                rawPen * cos(calcNormalizationR(IA_R, normalizationR));
            s.get_impact(i + j, impact::ePenDN) =
                rawPen * cos(calcNormalizationR(IAD_R, normalizationR));

            s.get_impact(i + j, impact::tToTarget) = tVec[j];
            s.get_impact(i + j, impact::tToTargetA) = tVec[j] / 3.1;
        }
    }

public:
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
        s.impactSize = (unsigned int)((max - min) / precision);
        s.impactSizeAligned = vSize - (s.impactSize % vSize) + s.impactSize;
        s.trajectories.resize(2 * s.impactSize);
        s.impactData.resize(impact::maxColumns * s.impactSizeAligned);

        if (nThreads > std::thread::hardware_concurrency()) {
            nThreads = std::thread::hardware_concurrency();
        }
        unsigned int length = ceil((double)s.impactSize / vSize);
        unsigned int assigned;
        if (length > nThreads) {
            assigned = nThreads;
        } else {
            assigned = length;
        }
        mtFunctionRunner(assigned, length, s.impactSize, this,
                         &shellCalc::multiTraj<AddTraj, Numerical, Hybrid>, &s);

        s.completedImpact = true;
    }

private:
    void checkRunImpact(shell &s) {
        if (!s.completedImpact) {
            std::cout << "Standard Not Calculated - Running automatically\n";
            calculateImpact<false, numerical::forwardEuler, false>(s);
        }
    }

    // Check Angles Section
    // template <short fusing> explanation
    // Possible Values:
    // 0 - Never Fusing
    // 1 - Check
    // 2 - Always Fusing
    enum fuseStatus { never, check, always };
    template <short fusing>
    void multiAngles(const unsigned int i, const double thickness,
                     const double inclination_R, const double fusingAngle,
                     shell *const shellPointer) {
        static_assert(fusing <= 2 && fusing >= 0, "Invalid fusing parameter");
        shell &s = *shellPointer;
        for (unsigned int j = 0; j < vSize; j++) {
            double fallAngleAdjusted =
                s.get_impact(i + j, impact::impactDataIndex::impactAHR) +
                inclination_R;
            double rawPen =
                s.get_impact(i + j, impact::impactDataIndex::rawPen);

            double penetrationCriticalAngle;
            penetrationCriticalAngle =
                (acos(thickness / rawPen) + s.get_normalizationR());
            penetrationCriticalAngle = std::isnan(penetrationCriticalAngle)
                                           ? 0
                                           : penetrationCriticalAngle;

            std::array<double, 4> criticalAngles = {
                s.get_ricochet0R(), s.get_ricochet1R(),
                penetrationCriticalAngle, fusingAngle};
            std::array<double, 4> out;
            for (unsigned int k = 0; k < 2; k++) {
                out[k] = acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                out[k] = std::isnan(out[k]) ? 0 : out[k];
            }

            // for (int k = 2; k < angle::maxColumns / 2; k++)
            {
                int k = angle::armor / 2;
                if (criticalAngles[k] < M_PI / 2) {
                    out[k] =
                        acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                    out[k] = std::isnan(out[k]) ? 0 : out[k];
                } else {
                    out[k] = M_PI / 2;
                }
            }
            {
                int k = angle::fuse / 2;
                if constexpr (fusing == fuseStatus::never) {
                    // std::cout << "0" << fusingAngle << "\n";
                    out[k] = M_PI / 2;
                } else if (fusing == fuseStatus::check) {
                    // std::cout << "1 " << fusingAngle << "\n";
                    out[k] =
                        acos(cos(criticalAngles[k]) / cos(fallAngleAdjusted));
                    out[k] = std::isnan(out[k]) ? 0 : out[k];
                } else if (fusing == fuseStatus::always) {
                    // std::cout << "2 " << fusingAngle << "\n";
                    out[k] = 0;
                }
            }

            s.get_angle(i + j, angle::angleDataIndex::ra0) = out[0];
            s.get_angle(i + j, angle::angleDataIndex::ra1) = out[1];
            s.get_angle(i + j, angle::angleDataIndex::armor) = out[2];
            s.get_angle(i + j, angle::angleDataIndex::fuse) = out[3];

            for (unsigned int k = 0; k < angle::maxColumns / 2; k++) {
                out[k] *= 180 / M_PI;
            }

            s.get_angle(i + j, angle::angleDataIndex::ra0D) = out[0];
            s.get_angle(i + j, angle::angleDataIndex::ra1D) = out[1];
            s.get_angle(i + j, angle::angleDataIndex::armorD) = out[2];
            s.get_angle(i + j, angle::angleDataIndex::fuseD) = out[3];
        }
    }

public:
    void calculateAngles(
        const double thickness, const double inclination, shell &s,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {

        checkRunImpact(s);

        s.angleData.resize(angle::maxColumns * s.impactSizeAligned);
        std::copy_n(s.get_impactPtr(0, impact::distance), s.impactSize,
                    s.get_anglePtr(0, angle::distance));

        unsigned int assigned;
        unsigned int length = (unsigned int)ceil((double)s.impactSize / vSize);
        if (length < nThreads) {
            assigned = length;
        } else {
            assigned = nThreads;
        }

        double inclination_R = inclination / 180 * M_PI;
        double fusingAngle;
        fusingAngle =
            acos(thickness / s.get_threshold()) + s.get_normalizationR();
        if (std::isnan(fusingAngle)) {
            mtFunctionRunner(assigned, length, s.impactSize, this,
                             &shellCalc::multiAngles<fuseStatus::always>,
                             thickness, inclination_R, fusingAngle, &s);
        } else {
            if (fusingAngle > M_PI / 2) {
                mtFunctionRunner(assigned, length, s.impactSize, this,
                                 &shellCalc::multiAngles<fuseStatus::never>,
                                 thickness, inclination_R, fusingAngle, &s);
            } else {
                mtFunctionRunner(assigned, length, s.impactSize, this,
                                 &shellCalc::multiAngles<fuseStatus::check>,
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
        if constexpr (fast) {
            double x = v_x * s.get_fuseTime();
            s.get_postPen(i, post::x, 0) = x;
            s.get_postPen(i, post::y, 0) = v_y * s.get_fuseTime();
            s.get_postPen(i, post::z, 0) = v_z * s.get_fuseTime();
            s.get_postPen(i, post::xwf, 0) =
                (thickness >= s.get_threshold()) * x +
                !(thickness >= s.get_threshold()) * -1;
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
                while (t < s.get_fuseTime()) {
                    for (int l = 0; l < 3; l++) {
                        pos[l] += velocities[l] * dtf;
                    }
                    // Calculate air density - likely unnecessary for this
                    // section as distances are so short
                    T = t0 - L * pos[1];
                    p = p0 * pow((1 - L * pos[1] / t0), (g * M / (R * L)));
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
                    dragIntermediary[0] = xz_dragIntermediary[0]; // x
                    dragIntermediary[1] =
                        (g - k * rho *
                                 (cw_1 * velocitiesSquared[1] +
                                  cw_2 * fabs(velocities[1])) *
                                 signum(velocities[1]));
                    dragIntermediary[2] = xz_dragIntermediary[1]; // z
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
                    (thickness >= s.get_threshold()) * pos[0] +
                    !(thickness >= s.get_threshold()) * -1;
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
            std::copy_n(s.get_impactPtr(distIndex, impact::impactAHR), vSize,
                        vAngleV);
            std::copy_n(s.get_impactPtr(distIndex, impact::rawPen), vSize,
                        penetrationV);
            std::copy_n(s.get_impactPtr(distIndex, impact::impactV), vSize,
                        v0V);
        } else {
            for (j = 0; (j + distIndex < s.impactSize) && (j < vSize); j++) {
                vAngleV[j] = s.get_impact(distIndex + j, impact::impactAHR);
                penetrationV[j] = s.get_impact(distIndex + j, impact::rawPen);
                v0V[j] = s.get_impact(distIndex + j, impact::impactV);
            }
            if (anglesIndex < s.postPenSize / s.impactSize) {
                for (; (j < vSize); j++) {
                    vAngleV[j] = s.get_impact(k, impact::impactAHR);
                    penetrationV[j] = s.get_impact(k, impact::rawPen);
                    v0V[j] = s.get_impact(k, impact::impactV);
                    k++;
                }
            }
        }

        for (unsigned int l = 0; l < vSize; l++) {
            double HA_R = hAngleV[l] * M_PI / 180;    // lateral  angle radians
            double VA_R = vAngleV[l] + inclination_R; // vertical angle radians
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
            std::copy_n(s->get_impactPtr(0, impact::distance), s->impactSize,
                        s->postPenData.begin() + s->postPenSize +
                            i * s->impactSize);
        }
        counter.fetch_add(1, std::memory_order_relaxed);
    }

    void parallelFillCopy(shell *const s, std::vector<double> *angles,
                          unsigned int nThreads) {
        std::vector<std::thread> threads;
        std::atomic<int> counter{0};
        unsigned int assigned;
        unsigned int length = angles->size();
        if (length < nThreads) {
            assigned = length;
        } else {
            assigned = nThreads;
        }
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
    void calculatePostPen(
        const double thickness, const double inclination, shell &s,
        std::vector<double> &angles, const bool changeDirection = true,
        const bool fast = false,
        const unsigned int nThreads = std::thread::hardware_concurrency()) {
        // Specifies whether normalization alters the trajectory of the shell
        // Though effect is not too significant either way
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
        /* Specifies whether to perform ballistic calculations or just multiply
           velocity by fusetime.
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
        unsigned int assigned;
        if (length < nThreads) {
            assigned = length;
        } else {
            assigned = nThreads;
        }
        mtFunctionRunner(assigned, length, s.postPenSize, this,
                         &shellCalc::multiPostPen<changeDirection, fast>,
                         thickness, inclination_R, &s);

        s.completedPostPen = true;
    }
};

} // namespace shell