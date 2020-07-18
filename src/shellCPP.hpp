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
    void multiTraj(const unsigned int start, shell &s,
                   double *vx, double *vy, double *tVec) {
        const double k = s.get_k();
        const double cw_2 = s.get_cw_2();
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

        auto delta = [&](double x, double &dx, 
                        double y, double &dy, 
                        double v_x, double &ddx, 
                        double v_y, double &ddy, bool update = false){
            update |= (y >= 0);
            double T, p, rho, dt_update = update * dt_min;
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
            //adamsBashforth5();
        }else{
            std::array<double, vSize> groupX;
            groupX.fill(0);
            std::array<double, vSize> groupY;
            for(unsigned int i=0, j=start; i<vSize; ++i, ++j){
                groupY[i] = j < s.impactSize ? 0: -1;
            }
            auto RK4Final = [](std::array<double, 4> &d) -> double{
                return (d[0] + 2*d[1] + 2*d[2] + d[3]) / 6;
            };
            auto rungeKutta4 = [&](std::size_t i){
                double &x = groupX[i], &y = groupY[i], 
                    &v_x = vx[i], &v_y = vy[i], &t = tVec[i];
                //double T, p, rho;
                bool update = (y >= 0); //Force update even if it becomes zero
                double dt_update = update * dt_min;
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
            };

            auto RK2Final = [](std::array<double, 2> &d) -> double{
                return (d[0] + d[1]) / 2;
            };
            auto rungeKutta2 = [&](std::size_t i){
                double &x = groupX[i], &y = groupY[i], 
                    &v_x = vx[i], &v_y = vy[i], &t = tVec[i];
                //double T, p, rho;
                double dt_update = (y >= 0) * dt_min;
                std::array<double, 2> dx, dy, ddx, ddy;
                
                delta(x, dx[0], y, dy[0], v_x, ddx[0], v_y, ddy[0]);
                delta(x+dx[0], dx[1], y+dy[0], dy[1], 
                    v_x+ddx[0], ddx[1], v_y+ddy[0], ddy[1], y >= 0); 
                //Force update even if it becomes zero

                x += RK2Final(dx); y += RK2Final(dy);
                v_x += RK2Final(ddx); v_y += RK2Final(ddy);
                t += dt_update;
            };

            auto forwardEuler = [&](std::size_t i){
                double &x = groupX[i], &y = groupY[i], 
                    &v_x = vx[i], &v_y = vy[i], &t = tVec[i];
                //double T, p, rho;
                double dt_update = (y >= 0) * dt_min;
                double dx, dy, ddx, ddy;

                delta(x, dx, y, dy, v_x, ddx, v_y, ddy);
                x += dx; y += dy;
                v_x += ddx; v_y += ddy;
                t += dt_update;
            };

            auto checkContinue = [&]() {
                bool any = false;
                for (unsigned int i = 0; i < vSize; ++i) {
                    any |= (groupY[i] >= 0);
                }
                return any;
            };
            while(checkContinue()){
                for(unsigned int i=0; i< vSize; ++i){
                    if constexpr(Numerical == numerical::forwardEuler){
                        forwardEuler(i);
                    }else if constexpr(Numerical == numerical::rungeKutta2){
                        rungeKutta2(i);
                    }else if constexpr(Numerical == numerical::rungeKutta4){
                        rungeKutta4(i);
                    }else{
                        static_assert(utility::falsy_v
                        <std::integral_constant<unsigned int, Numerical>>, 
                        "Invalid numerical algorithm");
                    }
                }
                if constexpr (AddTraj) {
                    for(unsigned int i=0, j=start; i<vSize; ++i, ++j){
                        if(j < s.impactSize){
                            s.trajectories[2 * (j)].push_back(groupX[i]);
                            s.trajectories[2 * (j) + 1].push_back(groupY[i]);
                        }
                    }
                }
            }
            for(std::size_t i=0; i< vSize; ++i){
                s.get_impact(start + i, impact::distance) = groupX[i];
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
        //for (unsigned int j = 0; (j + i < s.impactSize) & (j < vSize); j++) {
        //    singleTraj<AddTraj, Numerical, Hybrid>(i, j, s, vx, vy, tVec);
        //}
        multiTraj<AddTraj, Numerical>(i, s, vx, vy, tVec);
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