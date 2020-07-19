#ifndef _SHELL_WOWS_CALC_HPP_
#define _SHELL_WOWS_CALC_HPP_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <algorithm>
#include <array>
#include <atomic>

#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <type_traits>

#include "controlEnums.hpp"
#include "shell.hpp"
#include "utility.hpp"

namespace shell {
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
            }else{
                break;
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

    template <unsigned int Numerical>
    static constexpr bool isMultistep() {
        if constexpr (Numerical == numerical::adamsBashforth5) {
            return true;
        } else {
            return false;
        }
    }

    //https://godbolt.org/z/xj8qGz
    template <bool AddTraj, unsigned int Numerical>
    void multiTraj(const unsigned int &start, shell &s,
                   std::array<double, 2*vSize>& velocities, 
                   std::array<double, vSize>& tVec) {
        double k = s.get_k(), cw_2 = s.get_cw_2();
        if constexpr (AddTraj) {
            for(unsigned int i=0, j=start; i<vSize; ++i, ++j){
                if(j < s.impactSize){
                    s.trajectories[2 * (j)].clear();
                    s.trajectories[2 * (j) + 1].clear();
                    s.trajectories[2 * (j)].push_back(x0);
                    s.trajectories[2 * (j) + 1].push_back(y0);
                }
            }
        }

        std::array<double, 2*vSize> xy{};
        for(unsigned int i=0, j=start; i<vSize; ++i, ++j){
            xy[i+vSize] = j < s.impactSize ? 0: -1;
        }

        auto checkContinue = [&]() -> bool {
            bool any = false;
            for (unsigned int i = 0; i < vSize; ++i) {
                any |= (xy[i+vSize] >= 0);
            }
            return any;
        };

        auto delta = [&](const double &x, double &dx, 
                        double y, double &dy, 
                        const double &v_x, double &ddx, 
                        const double &v_y, double &ddy, 
                        bool update = false){
            update |= (y >= 0);
            double T, p, rho, dt_update = update ? dt_min:0;
            dx = dt_update * v_x;
            dy = dt_update * v_y;
            y += dy; //x not needed
            //Air Density
            T = t0 - L * y;
            p = p0 * pow(1 - L*y/ t0, gMRL);
            rho = p * M / (R * T);
            double kRho = k * rho;
            //Calculate Drag Components
            ddx = -1*dt_update*   kRho*(cw_1*v_x*v_x+cw_2*v_x);
            ddy = -1*dt_update*(g+kRho*(cw_1*v_y*v_y+cw_2*fabs(v_y)*signum(v_y)));
        };

        if constexpr (isMultistep<Numerical>()){
            if constexpr (Numerical == numerical::adamsBashforth5){
                std::array<double, 5*vSize> dx, dy, ddx, ddy;
                //0 -> vSize -> ... -> 5 * vSize
                unsigned int offset = 0; //Make it a circular buffer
                auto get = [&](const unsigned int &index, const unsigned int &stage) -> unsigned int{
                    return index + ((stage + offset) % 5) * vSize;
                };

                if(checkContinue()){ //AB1 == Euler Method - Length 1 Traj
                    auto ABF1 = [&](const std::array<double, 5*vSize> &d, const unsigned int &i, const bool &update){
                        return d[get(i, 0)] * update;
                    };
                    //ABG(0, ABF1);
                    for(unsigned int i=0; i<vSize; ++i){
                        double &x = xy[i], &y = xy[i+vSize], 
                        &v_x = velocities[i], &v_y = velocities[i+vSize], 
                        &t = tVec[i];
                        bool update = (y >= 0);
                        unsigned int index = get(i, 0);
                        delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y, ddy[index]);
                        x += ABF1(dx, i, update); y += ABF1(dy, i, update); 
                        v_x += ABF1(ddx, i, update); v_y += ABF1(ddy, i, update);
                        t += update * dt_min;
                    }
                    if constexpr (AddTraj){
                        for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                            s.trajectories[2 * (j)].push_back(xy[i]);
                            s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                        }
                    }
                    if(checkContinue()){ //2 AB2 - Length 2 Traj
                        auto ABF2 = [&](const std::array<double, 5*vSize> &d, const unsigned int &i, const bool &update){
                            return (3/2*d[get(i, 1)] - 1/2*d[get(i, 0)]) * update; 
                        };
                        for(unsigned int i=0; i<vSize; ++i){
                            double &x = xy[i], &y = xy[i+vSize], 
                            &v_x = velocities[i], &v_y = velocities[i+vSize], 
                            &t = tVec[i];
                            bool update = (y >= 0);
                            unsigned int index = get(i, 1);
                            delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y, ddy[index]);
                            x += ABF2(dx, i, update); y += ABF2(dy, i, update); 
                            v_x += ABF2(ddx, i, update); v_y += ABF2(ddy, i, update);
                            t += update * dt_min;
                        }
                        if constexpr (AddTraj){
                            for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                                s.trajectories[2 * (j)].push_back(xy[i]);
                                s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                            }
                        }
                        if(checkContinue()){ //3 AB3 - Length 3 Traj
                            auto ABF3 = [&](const std::array<double, 5*vSize> &d, const unsigned int &i, const bool &update){
                                return (23/12*d[get(i, 2)] - 16/12*d[get(i, 1)] + 5/12*d[get(i, 0)]) 
                                * update; 
                            };
                            for(unsigned int i=0; i<vSize; ++i){
                                double &x = xy[i], &y = xy[i+vSize], 
                                &v_x = velocities[i], &v_y = velocities[i+vSize], 
                                &t = tVec[i];
                                bool update = (y >= 0);
                                unsigned int index = get(i, 2);
                                delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y, ddy[index]);
                                x += ABF3(dx, i, update); y += ABF3(dy, i, update); 
                                v_x += ABF3(ddx, i, update); v_y += ABF3(ddy, i, update);
                                t += update * dt_min;
                            }
                            if constexpr (AddTraj){
                                for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                                    s.trajectories[2 * (j)].push_back(xy[i]);
                                    s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                                }
                            }
                            if(checkContinue()){ //4 AB4 - Length 4 Traj
                                auto ABF4 = [&](const std::array<double, 5*vSize> &d, const unsigned int &i, const bool &update){
                                    return (55/24*d[get(i, 3)] - 59/24*d[get(i, 2)] + 37/24*d[get(i, 1)]
                                    - 9/24*d[get(i, 0)]) * update; 
                                };
                                for(unsigned int i=0; i<vSize; ++i){
                                    double &x = xy[i], &y = xy[i+vSize], 
                                    &v_x = velocities[i], &v_y = velocities[i+vSize], 
                                    &t = tVec[i];
                                    bool update = (y >= 0);
                                    unsigned int index = get(i, 3);
                                    delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y, ddy[index]);
                                    x += ABF4(dx, i, update); y += ABF4(dy, i, update); 
                                    v_x += ABF4(ddx, i, update); v_y += ABF4(ddy, i, update);
                                    t += update * dt_min;
                                }
                                if constexpr (AddTraj){
                                    for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                                        s.trajectories[2 * (j)].push_back(xy[i]);
                                        s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                                    }
                                }
                                while(checkContinue()){ //5 AB5 - Length 5+ Traj
                                    auto ABF5 = [&](const std::array<double, 5*vSize> &d, const unsigned int &i, const bool &update){
                                        return (1901/720*d[get(i, 4)] - 2774/720*d[get(i, 3)] + 2616/720*d[get(i, 2)]
                                        - 1274/720*d[get(i, 1)] + 251/720*d[get(i, 0)]) * update; 
                                    };
                                    for(unsigned int i=0; i<vSize; ++i){
                                        double &x = xy[i], &y = xy[i+vSize], 
                                        &v_x = velocities[i], &v_y = velocities[i+vSize], 
                                        &t = tVec[i];
                                        bool update = (y >= 0);
                                        unsigned int index = get(i, 4);
                                        delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y, ddy[index]);
                                        x += ABF5(dx, i, update); y += ABF5(dy, i, update); 
                                        v_x += ABF5(ddx, i, update); v_y += ABF5(ddy, i, update);
                                        t += update * dt_min;
                                    }
                                    if constexpr (AddTraj){
                                        for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                                            s.trajectories[2 * (j)].push_back(xy[i]);
                                            s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                                        }
                                    }
                                    offset++; //Circle back 
                                    offset = offset == 5? 0 : offset;
                                }
                            }
                        }
                    }
                }
            }else{
                static_assert(utility::falsy_v
                <std::integral_constant<unsigned int, Numerical>>, 
                "Invalid multistep algorithm");
            }
        }else{
            auto RK4Final = [](std::array<double, 4> &d) -> double{
                return (std::fma(2, d[1], d[0]) + std::fma(2, d[2], d[3])) / 6;
            };

            auto RK2Final = [](std::array<double, 2> &d) -> double{
                return (d[0] + d[1]) / 2;
            };

            while(checkContinue()){
                for(unsigned int i=0; i< vSize; ++i){
                    double &x = xy[i], &y = xy[i+vSize],
                    &v_x = velocities[i], &v_y = velocities[i+vSize], &t = tVec[i];
                    if constexpr(Numerical == numerical::forwardEuler){
                        double dt_update = (y >= 0) ? dt_min:0;
                        double dx, dy, ddx, ddy;

                        delta(x, dx, y, dy, v_x, ddx, v_y, ddy);
                        x += dx; y += dy;
                        v_x += ddx; v_y += ddy;
                        t += dt_update;
                    }else if constexpr(Numerical == numerical::rungeKutta2){
                        double dt_update = (y >= 0) ? dt_min:0;
                        std::array<double, 2> dx, dy, ddx, ddy;
                        
                        delta(x, dx[0], y, dy[0], v_x, ddx[0], v_y, ddy[0]);
                        delta(x+dx[0], dx[1], y+dy[0], dy[1], 
                            v_x+ddx[0], ddx[1], v_y+ddy[0], ddy[1], y >= 0); 
                        //Force update even if it becomes zero

                        x += RK2Final(dx); y += RK2Final(dy);
                        v_x += RK2Final(ddx); v_y += RK2Final(ddy);
                        t += dt_update;
                    }else if constexpr(Numerical == numerical::rungeKutta4){
                        bool update = (y >= 0); //Force update even if it becomes zero
                        double dt_update = update ? dt_min:0;
                        std::array<double, 4> dx, dy, ddx, ddy;
                        // K1->K4
                        delta(x           , dx[0] , y           , dy[0], 
                            v_x         , ddx[0], v_y         , ddy[0]);
                        delta(x+dx[0]/2   , dx[1] , y+dy[0]/2   , dy[1], 
                            v_x+ddx[0]/2, ddx[1], v_y+ddy[0]/2, ddy[1], update); 
                        delta(x+dx[1]/2   , dx[2] , y+dy[1]/2   , dy[2], 
                            v_x+ddx[1]/2, ddx[2], v_y+ddy[1]/2, ddy[2], update);
                        delta(x+dx[2]     , dx[3] , y+dy[2]     , dy[3], 
                            v_x+ddx[2]  , ddx[3], v_y+ddy[2]  , ddy[3], update);

                        x += RK4Final(dx); y += RK4Final(dy);
                        v_x += RK4Final(ddx); v_y += RK4Final(ddy);
                        t += dt_update;
                    }else{
                        static_assert(utility::falsy_v
                        <std::integral_constant<unsigned int, Numerical>>, 
                        "Invalid single step algorithm");
                    }
                }
                if constexpr (AddTraj) {
                    for(unsigned int i=0, j=start; i<vSize & j < s.impactSize; ++i, ++j){
                        s.trajectories[2 * (j)].push_back(xy[i]);
                        s.trajectories[2 * (j) + 1].push_back(xy[i+vSize]);
                    }
                };
            }
        }  
        for(std::size_t i=0; i< vSize; ++i){
            s.get_impact(start + i, impact::distance) = xy[i];
        } 
    }

    // Several trajectories done in one chunk to allow for vectorization
    template <bool AddTraj, unsigned int Numerical, bool Hybrid, bool Fit,
              bool nonAP>
    void impactGroup(const unsigned int i, shell *const shellPointer) {
        shell &s = *shellPointer;
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_normalizationR();
        //std::cout<<"Entered\n";
        std::array<double, vSize * 2> velocities; 
        //0 -> (v_x) -> vSize -> (v_y) -> 2*vSize
        std::array<double, vSize> tVec{};
        for (unsigned int j = 0; j < vSize; j++) {
            if constexpr (!Fit) {
                s.get_impact(i + j, impact::launchAngle) =
                    precision * (i + j) + min;
            }
            double radianLaunch =
                s.get_impact(i + j, impact::launchAngle) * M_PI / 180;
            velocities[j] = s.get_v0() * cos(radianLaunch);
            velocities[j+vSize] = s.get_v0() * sin(radianLaunch);
        }
        //std::cout<<"Calculating\n";
        multiTraj<AddTraj, Numerical>(i, s, velocities, tVec);
        //std::cout<<"Processing\n";
        for (unsigned int j = 0; j < vSize; j++) {
            const double &v_x = velocities[j], &v_y = velocities[j+vSize];
            double IA_R = atan(v_y / v_x);

            s.get_impact(i + j, impact::impactAngleHorizontalRadians) = IA_R;
            double IAD_R = M_PI_2 + IA_R;
            double IA_D = IA_R * 180 / M_PI;
            s.get_impact(i + j, impact::impactAngleHorizontalDegrees) =
                IA_D * -1;
            s.get_impact(i + j, impact::impactAngleDeckDegrees) = 90 + IA_D;

            double IV = sqrt(v_x * v_x + v_y * v_y);
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
        //leave extra space to allow for copies into the region
        //ex: | 0, 1, 2, 3, 4, 5 | -> | 0, 1, 2, 3, 4, 5 |, 0, 1, 2 + padding
        //allows for easier vectorization of code that uses this data
        int processedSize = unalignedSize;
        if(processedSize % vSize != 0){
            processedSize += vSize - 1;
        }
    // Templates to reduce branching
        return vSize - (processedSize % vSize) + processedSize;
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
                         &shellCalc::impactGroup<false, Numerical, false, true>,
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
        if(thickness >= s.threshold){
            fusingAngle = 0;
        }else{
            fusingAngle = acos(thickness / s.threshold) + s.get_normalizationR();
        }
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

        std::copy_n(s.get_postPenPtr(i, post::angle, 0), 
                    std::min(vSize, s.postPenSize - i), hAngleV);
        std::copy_n(s.get_impactPtr(distIndex,
                                    impact::impactAngleHorizontalRadians),
                    vSize, vAngleV);
        std::copy_n(s.get_impactPtr(distIndex, impact::rawPenetration),
                    vSize, penetrationV);
        std::copy_n(s.get_impactPtr(distIndex, impact::impactVelocity),
                    vSize, v0V);

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
        //copies vSize - 1 from front - to avoid branching in hot loop
        std::copy_n(
            s.get_impactPtr(0, impact::impactAngleHorizontalRadians), vSize - 1, 
            s.get_impactPtr(s.impactSize, impact::impactAngleHorizontalRadians));
        std::copy_n(
            s.get_impactPtr(0, impact::impactVelocity), vSize - 1, 
            s.get_impactPtr(s.impactSize, impact::impactVelocity));
        std::copy_n(
            s.get_impactPtr(0, impact::rawPenetration), vSize - 1, 
            s.get_impactPtr(s.impactSize, impact::rawPenetration));
        /*std::cout<<vSize<<"\n";
        for(int i=0; i<s.impactSizeAligned; i++){
            for(int j=0; j<impact::maxColumns; j++){
                std::cout<<s.get_impact(i, j)<<" ";
            }
            std::cout<<"\n";
        }*/
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
#endif