#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "controlEnums.hpp"
#include "shell.hpp"
#include "utility.hpp"

namespace shell {
class shellCalc {
   private:
    // TODO: Static Constexpr these
    // Physical Constants     Description                  | Units
    double g = 9.8;        // Gravitational Constant       | m/(s^2)
    double t0 = 288.15;    // Temperature at Sea Level     | K
    double L = 0.0065;     // Atmospheric Lapse Rate       | C/m
    double p0 = 101325;    // Pressure at Sea Level        | Pa
    double R = 8.31447;    // Ideal Gas Constant           | J/(mol K)
    double M = 0.0289644;  // Molarity of Air at Sea Level | kg/mol
    double cw_1 = 1;

    double gMRL = (g * M) / (R * L);
    // Calculation Parameters
    double max = 25;        // Max Angle                    | degrees
    double min = 0;         // Min Angle                    | degrees
    double precision = .1;  // Angle Step                   | degrees
    double x0 = 0, y0 = 0;  // Starting x0, y0              | m
    double dt_min = .02;    // Time step                    | s

    static constexpr double timeMultiplier = 2.75;
    // For some reason the game has a different shell multiplier than the
    // global speed multiplier of 2.61 used for everything else.
    static constexpr double velocityPower = 1.4822064892953855;

    // delta t (dtf) for fusing needs to be smaller than the delta t (dt) used
    // for trajectories due to the shorter distances. Otherwise results become
    // jagged - precision suffers.
    double dtf = 0.0001;
    double xf0 = 0, yf0 = 0;

    // For vectorization - though probably not 100% necessary anymore since
    // intrinsics were removed [intrinsics had no significant improvements in
    // runtime] - but might still make it easier to vectorize
    static_assert(sizeof(double) == 8,
                  "Size of double is not 8 - required for vectorization");
    // Use float64 in the future
    static_assert(std::numeric_limits<double>::is_iec559,
                  "Type is not IEE754 compliant");
    // For setting up vectorization lane widths
    // 128bit because compilers cannot seem to generate 256bit instructions
    static constexpr uint32_t vSize = (128 / 8) / sizeof(double);
    static constexpr uint32_t minTasksPerThread = vSize;

   public:
    double calcNormalizationR(
        const double angle,
        const double normalizationR) {  // Input in radians
        return fabs(angle) > normalizationR ? fabs(angle) - normalizationR : 0;
        // Don't worry this branch goes away
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
    void mtFunctionRunnerSelected(uint32_t assigned, uint32_t length,
                                  uint32_t size, O object, F function,
                                  Args... args) {
        if constexpr (multiThreaded) {
            std::atomic<int> counter{0};
            // moodycamel::ConcurrentQueue<int> workQueue;
            std::vector<std::thread> threads(assigned - 1);
            for (uint32_t i = 0; i < assigned - 1; i++) {
                threads[i] = std::thread([&, i] {
                    mtWorker(counter, length, i, object, function, args...);
                });
            }

            mtWorker(counter, length, assigned - 1, object, function, args...);

            for (uint32_t i = 0; i < assigned - 1; i++) {
                threads[i].join();
            }
        } else {
            for (uint32_t i = 0; i < size; i += vSize) {
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
            if (index < length) {
                // std::cout<<index<<"\n";
                (object->*function)(index * vSize, args...);
            }
        }
    }

    uint32_t assignThreadNum(uint32_t length, uint32_t nThreads) {
        if (length > nThreads * minTasksPerThread) {
            return nThreads;
        } else {
            return ceil(static_cast<double>(length) / minTasksPerThread);
        }
    }

    // https://godbolt.org/z/4b1sn5
    template <bool AddTraj, numerical Numerical>
    void multiTraj(const uint32_t &start, shell &s,
                   std::array<double, 3 * vSize> &velocities) {
        const double k = s.get_k(), cw_2 = s.get_cw_2();
        if constexpr (AddTraj) {
            for (uint32_t i = 0, j = start; i < vSize; ++i, ++j) {
                if (j < s.impactSize) {
                    s.trajectories[2 * (j)].clear();
                    s.trajectories[2 * (j) + 1].clear();
                    s.trajectories[2 * (j)].push_back(x0);
                    s.trajectories[2 * (j) + 1].push_back(y0);
                }
            }
        }

        std::array<double, 2 * vSize> xy{};
        for (uint32_t i = 0, j = start; i < vSize; ++i, ++j) {
            xy[i + vSize] = j < s.impactSize ? 0 : -1;
        }

        auto checkContinue = [&]() -> bool {
            bool any = false;
            for (uint32_t i = 0; i < vSize; ++i) {
                any |= (xy[i + vSize] >= 0);
            }
            return any;
        };

        // Helpers
        auto delta = [&](const double x, double &dx, double y, double &dy,
                         const double v_x, double &ddx, const double v_y,
                         double &ddy, bool update = false) {
            update |= (y >= 0);
            double T, p, rho, dt_update = update * dt_min;
            dx = dt_update * v_x;
            dy = dt_update * v_y;
            y += dy;  // x not needed
            // Air Density
            T = t0 - L * y;
            p = p0 * pow(1 - L * y / t0, gMRL);
            rho = p * M / (R * T);
            double kRho = k * rho;
            // Calculate Drag Components
            double velocityMagnitude = sqrt(v_x * v_x + v_y * v_y);
            ddx = -1 * dt_update * kRho *
                  (cw_1 * v_x * velocityMagnitude + cw_2 * v_x);
            ddy = -1 * dt_update *
                  (g + kRho * (cw_1 * v_y * velocityMagnitude +
                               cw_2 * fabs(v_y) * signum(v_y)));
        };

        auto getIntermediate = [](uint32_t index, uint32_t stage) {
            return index + stage * vSize;
        };

        auto RK4Final = [&](std::array<double, 4 * vSize> &d,
                            uint32_t index) -> double {
            // Adds deltas in Runge Kutta 4 manner
            return (std::fma(2, d[getIntermediate(index, 1)],
                             d[getIntermediate(index, 0)]) +
                    std::fma(2, d[getIntermediate(index, 2)],
                             d[getIntermediate(index, 3)])) /
                   6;
        };

        auto RK2Final = [&](std::array<double, 2 * vSize> &d,
                            uint32_t index) -> double {
            // Adds deltas in Runge Kutta 2 manner
            return (d[getIntermediate(index, 0)] +
                    d[getIntermediate(index, 1)]) /
                   2;
        };
        // TODO: Add numerical orders
        if constexpr (isMultistep<Numerical>()) {
            if constexpr (Numerical == numerical::adamsBashforth5) {
                std::array<double, 5 * vSize> dx, dy, ddx, ddy;
                // 0 -> vSize -> ... -> 5 * vSize
                uint32_t offset = 0;  // Make it a circular buffer
                auto get = [&](const uint32_t &index,
                               const uint32_t &stage) -> uint32_t {
                    return index + ((stage + offset) % 5) * vSize;
                };

                // Fill in first 5 w/ RK2
                std::array<double, 2 * vSize> rdx, rdy, rddx, rddy;
                for (int stage = 0; (stage < 4) & checkContinue(); ++stage) {
                    for (uint32_t i = 0; i < vSize; ++i) {
                        double &x = xy[i], &y = xy[i + vSize],
                               &v_x = velocities[i],
                               &v_y = velocities[i + vSize],
                               &t = velocities[i + vSize * 2];

                        // RK2
                        double dt_update = (y >= 0) * dt_min;
                        // std::array<double, 2> rdx, rdy, rddx, rddy;

                        delta(x, rdx[getIntermediate(i, 0)], y,
                              rdy[getIntermediate(i, 0)], v_x,
                              rddx[getIntermediate(i, 0)], v_y,
                              rddy[getIntermediate(i, 0)]);
                        delta(x + rdx[getIntermediate(i, 0)],
                              rdx[getIntermediate(i, 1)],
                              y + rdy[getIntermediate(i, 0)],
                              rdy[getIntermediate(i, 1)],
                              v_x + rddx[getIntermediate(i, 0)],
                              rddx[getIntermediate(i, 1)],
                              v_y + rddy[getIntermediate(i, 0)],
                              rddy[getIntermediate(i, 1)], y >= 0);
                        // Force update even if it becomes zero

                        double fdx = RK2Final(rdx, i), fdy = RK2Final(rdy, i),
                               fddx = RK2Final(rddx, i),
                               fddy = RK2Final(rddy, i);
                        x += fdx;
                        y += fdy;
                        v_x += fddx;
                        v_y += fddy;
                        dx[get(i, stage)] = fdx;
                        dy[get(i, stage)] = fdy;
                        ddx[get(i, stage)] = fddx;
                        ddy[get(i, stage)] = fddy;
                        t += dt_update;
                    }
                    if constexpr (AddTraj) {
                        const uint32_t loopSize =
                            std::min<uint32_t>(vSize, s.impactSize - start);
                        for (uint32_t i = 0, j = start; i < loopSize;
                             ++i, ++j) {
                            s.trajectories[2 * (j)].push_back(xy[i]);
                            s.trajectories[2 * (j) + 1].push_back(
                                xy[i + vSize]);
                        }
                    }
                }

                while (checkContinue()) {  // 5 AB5 - Length 5+ Traj
                    auto ABF5 = [&](const std::array<double, 5 * vSize> &d,
                                    const uint32_t &i, const bool &update) {
                        // Adds deltas in Adams Bashforth 5 manner
                        return (1901 / 720 * d[get(i, 4)] -
                                2774 / 720 * d[get(i, 3)] +
                                2616 / 720 * d[get(i, 2)] -
                                1274 / 720 * d[get(i, 1)] +
                                251 / 720 * d[get(i, 0)]) *
                               update;
                    };
                    for (uint32_t i = 0; i < vSize; ++i) {
                        double &x = xy[i], &y = xy[i + vSize],
                               &v_x = velocities[i],
                               &v_y = velocities[i + vSize],
                               &t = velocities[i + vSize * 2];
                        bool update = (y >= 0);
                        uint32_t index = get(i, 4);  // Write to index
                        delta(x, dx[index], y, dy[index], v_x, ddx[index], v_y,
                              ddy[index]);
                        x += ABF5(dx, i, update);
                        y += ABF5(dy, i, update);
                        v_x += ABF5(ddx, i, update);
                        v_y += ABF5(ddy, i, update);
                        t += update * dt_min;
                    }
                    if constexpr (AddTraj) {
                        const uint32_t loopSize =
                            std::min<uint32_t>(vSize, s.impactSize - start);
                        for (uint32_t i = 0, j = start; i < loopSize;
                             ++i, ++j) {
                            s.trajectories[2 * (j)].push_back(xy[i]);
                            s.trajectories[2 * (j) + 1].push_back(
                                xy[i + vSize]);
                        }
                    }
                    offset++;  // Circle back
                    offset %= 5;
                }
            } else {
                static_assert(utility::falsy_v<std::integral_constant<
                                  uint32_t, toUnderlying(Numerical)>>,
                              "Invalid multistep algorithm");
            }
        } else {
            while (checkContinue()) {
                if constexpr (Numerical == numerical::forwardEuler) {
                    std::array<double, vSize> dx, dy, ddx, ddy;
                    for (uint32_t i = 0; i < vSize; ++i) {
                        double &x = xy[i], &y = xy[i + vSize],
                               &v_x = velocities[i],
                               &v_y = velocities[i + vSize],
                               &t = velocities[i + vSize * 2];
                        double dt_update = (y >= 0) * dt_min;
                        // double dx, dy, ddx, ddy;

                        delta(x, dx[i], y, dy[i], v_x, ddx[i], v_y, ddy[i]);
                        x += dx[i];
                        y += dy[i];
                        v_x += ddx[i];
                        v_y += ddy[i];
                        t += dt_update;
                    }
                } else if constexpr (Numerical == numerical::rungeKutta2) {
                    std::array<double, 2 * vSize> dx, dy, ddx, ddy;
                    for (uint32_t i = 0; i < vSize; ++i) {
                        double &x = xy[i], &y = xy[i + vSize],
                               &v_x = velocities[i],
                               &v_y = velocities[i + vSize],
                               &t = velocities[i + vSize * 2];
                        double dt_update = (y >= 0) * dt_min;
                        // std::array<double, 2> dx, dy, ddx, ddy;

                        delta(x, dx[getIntermediate(i, 0)], y,
                              dy[getIntermediate(i, 0)], v_x,
                              ddx[getIntermediate(i, 0)], v_y,
                              ddy[getIntermediate(i, 0)]);
                        delta(x + dx[getIntermediate(i, 0)],
                              dx[getIntermediate(i, 1)],
                              y + dy[getIntermediate(i, 0)],
                              dy[getIntermediate(i, 1)],
                              v_x + ddx[getIntermediate(i, 0)],
                              ddx[getIntermediate(i, 1)],
                              v_y + ddy[getIntermediate(i, 0)],
                              ddy[getIntermediate(i, 1)], y >= 0);
                        // Force update even if it becomes zero

                        x += RK2Final(dx, i);
                        y += RK2Final(dy, i);
                        v_x += RK2Final(ddx, i);
                        v_y += RK2Final(ddy, i);
                        t += dt_update;
                    }
                } else if constexpr (Numerical == numerical::rungeKutta4) {
                    std::array<double, 4 * vSize> dx, dy, ddx, ddy;
                    for (uint32_t i = 0; i < vSize; ++i) {
                        double &x = xy[i], &y = xy[i + vSize],
                               &v_x = velocities[i],
                               &v_y = velocities[i + vSize],
                               &t = velocities[i + vSize * 2];
                        bool update =
                            (y >= 0);  // Force update even if it becomes zero
                        double dt_update = update * dt_min;
                        // std::array<double, 4> dx, dy, ddx, ddy;
                        // K1->K4
                        delta(x, dx[getIntermediate(i, 0)], y,
                              dy[getIntermediate(i, 0)], v_x,
                              ddx[getIntermediate(i, 0)], v_y,
                              ddy[getIntermediate(i, 0)]);
                        delta(x + dx[getIntermediate(i, 0)] / 2,
                              dx[getIntermediate(i, 1)],
                              y + dy[getIntermediate(i, 0)] / 2,
                              dy[getIntermediate(i, 1)],
                              v_x + ddx[getIntermediate(i, 0)] / 2,
                              ddx[getIntermediate(i, 1)],
                              v_y + ddy[getIntermediate(i, 0)] / 2,
                              ddy[getIntermediate(i, 1)], update);
                        delta(x + dx[getIntermediate(i, 1)] / 2,
                              dx[getIntermediate(i, 2)],
                              y + dy[getIntermediate(i, 1)] / 2,
                              dy[getIntermediate(i, 2)],
                              v_x + ddx[getIntermediate(i, 1)] / 2,
                              ddx[getIntermediate(i, 2)],
                              v_y + ddy[getIntermediate(i, 1)] / 2,
                              ddy[getIntermediate(i, 2)], update);
                        delta(x + dx[getIntermediate(i, 2)],
                              dx[getIntermediate(i, 3)],
                              y + dy[getIntermediate(i, 2)],
                              dy[getIntermediate(i, 3)],
                              v_x + ddx[getIntermediate(i, 2)],
                              ddx[getIntermediate(i, 3)],
                              v_y + ddy[getIntermediate(i, 2)],
                              ddy[getIntermediate(i, 3)], update);

                        x += RK4Final(dx, i);
                        y += RK4Final(dy, i);
                        v_x += RK4Final(ddx, i);
                        v_y += RK4Final(ddy, i);
                        t += dt_update;
                    }
                } else {
                    static_assert(utility::falsy_v<std::integral_constant<
                                      uint32_t, toUnderlying(Numerical)>>,
                                  "Invalid single step algorithm");
                }

                if constexpr (AddTraj) {
                    const uint32_t loopSize =
                        std::min<uint32_t>(vSize, s.impactSize - start);
                    for (uint32_t i = 0, j = start; i < loopSize;
                         ++i, ++j) {  // Breaks Vectorization
                        s.trajectories[2 * (j)].push_back(xy[i]);
                        s.trajectories[2 * (j) + 1].push_back(xy[i + vSize]);
                    }
                };
            }
        }
        std::copy_n(xy.begin(), vSize,
                    s.get_impactPtr(start, impact::impactIndices::distance));
    }

    // Several trajectories done in one chunk to allow for vectorization
    template <bool AddTraj, numerical Numerical, bool Hybrid, bool Fit,
              bool nonAP>
    void impactGroup(const uint32_t i, shell &s) {
        // shell &s = *shellPointer;
        const double pPPC = s.get_pPPC();
        const double normalizationR = s.get_normalizationR();
        // std::cout<<"Entered\n";
        std::array<double, vSize * 3> velocitiesTime{};
        // 0 -> (v_x) -> vSize -> (v_y) -> 2*vSize -> (t) -> 3*vSize
        for (uint32_t j = 0; j < vSize; j++) {
            double radianLaunch;
            if constexpr (!Fit) {
                double degreeLaunch = precision * (i + j) + min;
                s.get_impact(i + j, impact::impactIndices::launchAngle) =
                    degreeLaunch;
                radianLaunch = degreeLaunch * M_PI / 180;
            } else {
                radianLaunch =
                    s.get_impact(i + j, impact::impactIndices::launchAngle) *
                    M_PI / 180;
            }
            velocitiesTime[j] = s.get_v0() * cos(radianLaunch);
            velocitiesTime[j + vSize] = s.get_v0() * sin(radianLaunch);
        }
        // std::cout<<"Calculating\n";
        multiTraj<AddTraj, Numerical>(i, s, velocitiesTime);
        // std::cout<<"Processing\n";
        for (uint32_t j = 0; j < vSize; j++) {
            const double &v_x = velocitiesTime[j],
                         &v_y = velocitiesTime[j + vSize];
            double IA_R = atan(v_y / v_x);

            s.get_impact(i + j,
                         impact::impactIndices::impactAngleHorizontalRadians) =
                IA_R;
            double IAD_R = M_PI_2 + IA_R;
            double IA_D = IA_R * 180 / M_PI;
            s.get_impact(i + j,
                         impact::impactIndices::impactAngleHorizontalDegrees) =
                IA_D * -1;
            s.get_impact(i + j, impact::impactIndices::impactAngleDeckDegrees) =
                90 + IA_D;

            double IV = sqrt(v_x * v_x + v_y * v_y);
            s.get_impact(i + j, impact::impactIndices::impactVelocity) = IV;

            double time = velocitiesTime[j + 2 * vSize];
            s.get_impact(i + j, impact::impactIndices::timeToTarget) = time;
            s.get_impact(i + j, impact::impactIndices::timeToTargetAdjusted) =
                time / timeMultiplier;

            if constexpr (!Fit) {
                if constexpr (nonAP) {
                    s.get_impact(i + j, impact::impactIndices::rawPenetration) =
                        s.nonAP;
                    s.get_impact(
                        i + j,
                        impact::impactIndices::effectivePenetrationHorizontal) =
                        s.nonAP;
                    s.get_impact(
                        i + j,
                        impact::impactIndices::effectivePenetrationDeck) =
                        s.nonAP;
                    s.get_impact(i + j,
                                 impact::impactIndices::
                                     effectivePenetrationHorizontalNormalized) =
                        s.nonAP;
                    s.get_impact(i + j,
                                 impact::impactIndices::
                                     effectivePenetrationDeckNormalized) =
                        s.nonAP;
                } else {
                    // double rawPenetration = pPPC * pow(IV, 1.1001);
                    // double rawPenetration = pPPC * pow(IV, 1.38);
                    // double rawPenetration = pPPC * pow(IV, 1.53803192);
                    // double rawPenetration = pPPC * pow(IV, 1.54562941);
                    double rawPenetration = pPPC * pow(IV, velocityPower);
                    s.get_impact(i + j, impact::impactIndices::rawPenetration) =
                        rawPenetration;

                    s.get_impact(
                        i + j,
                        impact::impactIndices::effectivePenetrationHorizontal) =
                        rawPenetration * cos(IA_R);
                    s.get_impact(
                        i + j,
                        impact::impactIndices::effectivePenetrationDeck) =
                        rawPenetration * cos(IAD_R);

                    s.get_impact(i + j,
                                 impact::impactIndices::
                                     effectivePenetrationHorizontalNormalized) =
                        rawPenetration *
                        cos(calcNormalizationR(IA_R, normalizationR));
                    s.get_impact(i + j,
                                 impact::impactIndices::
                                     effectivePenetrationDeckNormalized) =
                        rawPenetration *
                        cos(calcNormalizationR(IAD_R, normalizationR));
                }
            }
        }
    }

   public:
    uint32_t calculateAlignmentSize(uint32_t unalignedSize) {
        // leave extra space to allow for copies into the region
        // ex: | 0, 1, 2, 3, 4, 5 | -> | 0, 1, 2, 3, 4, 5 |, 0, 1, 2 + padding
        // allows for easier vectorization of code that uses this data
        int processedSize = unalignedSize;
        if (processedSize % vSize != 0) {
            processedSize += vSize - 1;
        }
        // Templates to reduce branching
        return vSize - (processedSize % vSize) + processedSize;
    }
    // Templates to reduce branching
    template <auto Numerical, bool Hybrid>
    void calculateImpact(
        shell &s, bool addTraj,
        uint32_t nThreads = std::thread::hardware_concurrency()) {
        if (addTraj) {
            calculateImpact<true, Numerical, Hybrid>(s, nThreads);
        } else {
            calculateImpact<false, Numerical, Hybrid>(s, nThreads);
        }
    }

    template <bool AddTraj, auto Numerical, bool Hybrid>
    void calculateImpact(
        shell &s, uint32_t nThreads = std::thread::hardware_concurrency()) {
        if (s.enableNonAP) {
            calculateImpact<true, Numerical, Hybrid, true>(s, nThreads);
        } else {
            calculateImpact<true, Numerical, Hybrid, false>(s, nThreads);
        }
    }

    template <bool AddTraj, auto Numerical, bool Hybrid, bool nonAP>
    void calculateImpact(
        shell &s, uint32_t nThreads = std::thread::hardware_concurrency()) {
        s.impactSize = ((max / precision - min / precision)) + 1;
        s.impactSizeAligned = calculateAlignmentSize(s.impactSize);
        if constexpr (AddTraj) {
            s.trajectories.resize(2 * s.impactSize);
        }
        s.impactData.resize(impact::maxColumns * s.impactSizeAligned);

        if (nThreads > std::thread::hardware_concurrency()) {
            nThreads = std::thread::hardware_concurrency();
        }
        uint32_t length = ceil(static_cast<double>(s.impactSize) / vSize);
        uint32_t assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(
            assigned, length, s.impactSize, this,
            &shellCalc::impactGroup<AddTraj, Numerical, Hybrid, false, nonAP>,
            std::ref(s));

        s.completedImpact = true;
    }

    template <auto Numerical>
    void calculateFit(shell &s,
                      uint32_t nThreads = std::thread::hardware_concurrency()) {
        if (nThreads > std::thread::hardware_concurrency()) {
            nThreads = std::thread::hardware_concurrency();
        }
        uint32_t length = ceil(static_cast<double>(s.impactSize) / vSize);
        uint32_t assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(
            assigned, length, s.impactSize, this,
            &shellCalc::impactGroup<false, Numerical, false, true, false>, &s);
    }

   private:
    void checkRunImpact(shell &s) {
        if (!s.completedImpact) {
            std::cout << "Standard Not Calculated - Running automatically\n";
            calculateImpact<false, numerical::forwardEuler, false>(s);
        }
    }

    // Check Angles Section
    // template <fuseStatus fusing> explanation:
    // Fusing is done using templates to reduce in loop branching and
    // computational time in some cases.

    // Possible Values: 0 - Never Fusing 1 - Check 2 - Always Fusing
    enum class fuseStatus { never, check, always };
    template <fuseStatus fusing, bool nonAP, bool nonAPPerforated,
              bool disableRicochet>
    void multiAngles(const uint32_t i, const double thickness,
                     const double inclination_R, const double fusingAngle,
                     shell &s) {
        static_assert(toUnderlying(fusing) <= 2 && toUnderlying(fusing) >= 0,
                      "Invalid fusing parameter");
        // shell &s = *shellPointer;
        const uint32_t ISA = s.impactSizeAligned;
        // ^^^
        // Otherwise Clang would think that assigned values are
        // "value[s] that could not be identified as reduction is used outside
        // the loop" This doesn't vectorize anyways - because of the acos's -
        // but that's there so that when acos vectorization is added
        // to compilers this will autovectorize
        for (uint32_t j = 0; j < vSize; j++) {
            double fallAngleAdjusted =
                s.impactData[i + j +
                             toUnderlying(impact::impactIndices::
                                              impactAngleHorizontalRadians) *
                                 ISA] +
                inclination_R;
            double rawPenetration =
                s.impactData[i + j +
                             toUnderlying(
                                 impact::impactIndices::rawPenetration) *
                                 ISA];

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

                penetrationCriticalAngle =
                    thickness > rawPenetration ? 0 : penetrationCriticalAngle;
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
            for (uint32_t k = 0; k < 2; k++) {
                double quotient =
                    cos(criticalAngles[k]) / cos(fallAngleAdjusted);
                out[k] = acos(quotient);
                out[k] = fabs(quotient) > 1 ? 0 : out[k];
            }

            {
                int k = toUnderlying(angle::angleIndices::armorRadians) / 2;
                double quotient =
                    cos(criticalAngles[k]) / cos(fallAngleAdjusted);
                out[k] = acos(quotient);
                out[k] = fabs(quotient) > 1 ? 0 : out[k];
                out[k] = criticalAngles[k] < M_PI_2 ? out[k] : M_PI_2;
                // Can't use ifs because for some reason (cond) ? (v1) : (v2)
                // is not equal to if(cond) v1 else v2 - creates jumps
            }
            {
                int k = toUnderlying(angle::angleIndices::fuseRadians) / 2;
                if constexpr (fusing == fuseStatus::never) {
                    out[k] = M_PI_2;
                } else if (fusing == fuseStatus::check) {
                    double quotient =
                        cos(criticalAngles[k]) / cos(fallAngleAdjusted);
                    out[k] = acos(quotient);
                    out[k] = fabs(quotient) > 1 ? 0 : out[k];
                } else if (fusing == fuseStatus::always) {
                    out[k] = 0;
                }
            }

            s.get_angle(i + j, angle::angleIndices::ricochetAngle0Radians) =
                out[0];
            s.get_angle(i + j, angle::angleIndices::ricochetAngle1Radians) =
                out[1];
            s.get_angle(i + j, angle::angleIndices::armorRadians) = out[2];
            s.get_angle(i + j, angle::angleIndices::fuseRadians) = out[3];

            for (uint32_t k = 0; k < angle::maxColumns / 2; k++) {
                out[k] *= 180 / M_PI;
            }

            s.get_angle(i + j, angle::angleIndices::ricochetAngle0Degrees) =
                out[0];
            s.get_angle(i + j, angle::angleIndices::ricochetAngle1Degrees) =
                out[1];
            s.get_angle(i + j, angle::angleIndices::armorDegrees) = out[2];
            s.get_angle(i + j, angle::angleIndices::fuseDegrees) = out[3];
        }
    }

   public:
    void calculateAngles(
        const double thickness, const double inclination, shell &s,
        const uint32_t nThreads = std::thread::hardware_concurrency()) {
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
        const uint32_t nThreads = std::thread::hardware_concurrency()) {
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
        const uint32_t nThreads = std::thread::hardware_concurrency()) {
        checkRunImpact(s);

        s.angleData.resize(angle::maxColumns * s.impactSizeAligned);
        std::copy_n(s.get_impactPtr(0, impact::impactIndices::distance),
                    s.impactSize,
                    s.get_anglePtr(0, angle::angleIndices::distance));

        uint32_t length = static_cast<uint32_t>(
            ceil(static_cast<double>(s.impactSize) / vSize));
        uint32_t assigned = assignThreadNum(length, nThreads);

        double inclination_R = inclination / 180 * M_PI;
        double fusingAngle;
        if (thickness >= s.threshold) {
            fusingAngle = 0;
        } else {
            fusingAngle =
                acos(thickness / s.threshold) + s.get_normalizationR();
        }
        if (thickness > s.threshold) {
            mtFunctionRunner(
                assigned, length, s.impactSize, this,
                &shellCalc::multiAngles<fuseStatus::always, nonAP,
                                        nonAPPerforated, disableRicochet>,
                thickness, inclination_R, fusingAngle, std::ref(s));
        } else {
            if (fusingAngle > M_PI_2) {
                mtFunctionRunner(
                    assigned, length, s.impactSize, this,
                    &shellCalc::multiAngles<fuseStatus::never, nonAP,
                                            nonAPPerforated, disableRicochet>,
                    thickness, inclination_R, fusingAngle, std::ref(s));
            } else {
                mtFunctionRunner(
                    assigned, length, s.impactSize, this,
                    &shellCalc::multiAngles<fuseStatus::check, nonAP,
                                            nonAPPerforated, disableRicochet>,
                    thickness, inclination_R, fusingAngle, std::ref(s));
            }
        }
        s.completedAngles = true;
    }

    // Post-Penetration Section

   private:
    template <bool fast>
    void postPenTraj(const uint32_t i, shell &s, double v_x, double v_y,
                     double v_z, double thickness) {
        const double notFusedCode = -1;
        if constexpr (fast) {
            bool positive = v_x > 0;
            double x = v_x * s.fuseTime * positive;
            s.get_postPen(i, post::postPenIndices::x, 0) = x;
            s.get_postPen(i, post::postPenIndices::y, 0) =
                v_y * s.fuseTime * positive;
            s.get_postPen(i, post::postPenIndices::z, 0) =
                v_z * s.fuseTime * positive;

            bool fuse = thickness >= s.threshold;
            s.get_postPen(i, post::postPenIndices::xwf, 0) =
                (fuse)*x + !(fuse)*notFusedCode;
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
                s.get_postPen(i, post::postPenIndices::x, 0) = pos[0];
                s.get_postPen(i, post::postPenIndices::y, 0) = pos[1];
                s.get_postPen(i, post::postPenIndices::z, 0) = pos[2];
                s.get_postPen(i, post::postPenIndices::xwf, 0) =
                    (thickness >= s.threshold) * pos[0] +
                    !(thickness >= s.threshold) * notFusedCode;
            } else {
                s.get_postPen(i, post::postPenIndices::x, 0) = 0;
                s.get_postPen(i, post::postPenIndices::y, 0) = 0;
                s.get_postPen(i, post::postPenIndices::z, 0) = 0;
                s.get_postPen(i, post::postPenIndices::xwf, 0) = 0;
            }
        }
    }

    template <bool changeDirection, bool fast>
    void multiPostPen(uint32_t i, const double thickness,
                      const double inclination_R, shell &s) {
        // shell &s = *shellPointer;
        double hAngleV[vSize], vAngleV[vSize];
        double v0V[vSize], penetrationV[vSize], eThicknessV[vSize];
        double v_x[vSize], v_y[vSize], v_z[vSize];
        uint32_t distIndex = (i < s.impactSize) ? i : i % s.impactSize;

        std::copy_n(s.get_postPenPtr(i, post::postPenIndices::angle, 0),
                    std::min<uint32_t>(vSize, s.postPenSize - i), hAngleV);
        std::copy_n(
            s.get_impactPtr(
                distIndex, impact::impactIndices::impactAngleHorizontalRadians),
            vSize, vAngleV);
        std::copy_n(
            s.get_impactPtr(distIndex, impact::impactIndices::rawPenetration),
            vSize, penetrationV);
        std::copy_n(
            s.get_impactPtr(distIndex, impact::impactIndices::impactVelocity),
            vSize, v0V);

        for (uint32_t l = 0; l < vSize; l++) {
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

        //#pragma clang loop vectorize(enable)
        const uint32_t loopSize = std::min<uint32_t>(vSize, s.postPenSize - i);
        for (uint32_t j = 0; j < loopSize; j++) {
            postPenTraj<fast>(i + j, s, v_x[j], v_y[j], v_z[j], eThicknessV[j]);
        }
        // std::cout<<index<<" Completed\n";
    }

    // Probably unnecessary...
    void fillCopy(std::atomic<int> &counter, uint32_t assigned, uint32_t id,
                  shell *const s, std::vector<double> *angles) {
        for (uint32_t i = angles->size() * id / assigned;
             i < angles->size() * (id + 1) / assigned; i++) {
            std::fill_n(s->get_postPenPtr(0, post::postPenIndices::angle, i),
                        s->impactSize, static_cast<double>((*angles)[i]));
            std::copy_n(
                s->get_impactPtr(0, impact::impactIndices::distance),
                s->impactSize,
                s->postPenData.begin() + s->postPenSize + i * s->impactSize);
        }
        counter.fetch_add(1, std::memory_order_relaxed);
    }

    void parallelFillCopy(shell *const s, std::vector<double> *angles,
                          uint32_t nThreads) {
        std::vector<std::thread> threads;
        std::atomic<int> counter{0};
        uint32_t length = angles->size();
        uint32_t assigned = assignThreadNum(length, nThreads);
        for (uint32_t i = 0; i < assigned - 1; i++) {
            threads.emplace_back(
                [&, i] { fillCopy(counter, assigned, i, s, angles); });
        }
        fillCopy(counter, assigned, assigned - 1, s, angles);
        for (uint32_t i = 0; i < assigned - 1; i++) {
            threads[i].join();
        }
    }

   public:
    // Again templates to reduce branching
    void calculatePostPen(
        const double thickness, const double inclination, shell &s,
        std::vector<double> &angles, const bool changeDirection = true,
        const bool fast = false,
        const uint32_t nThreads = std::thread::hardware_concurrency()) {
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
        const uint32_t nThreads = std::thread::hardware_concurrency()) {
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
                          const uint32_t nThreads) {
        checkRunImpact(s);

        s.postPenSize = s.impactSize * angles.size();
        s.postPenData.resize(6 * s.postPenSize);

        parallelFillCopy(&s, &angles, nThreads);
        // copies vSize - 1 from front - to avoid branching in hot loop
        std::copy_n(s.get_impactPtr(
                        0, impact::impactIndices::impactAngleHorizontalRadians),
                    vSize - 1,
                    s.get_impactPtr(
                        s.impactSize,
                        impact::impactIndices::impactAngleHorizontalRadians));
        std::copy_n(s.get_impactPtr(0, impact::impactIndices::impactVelocity),
                    vSize - 1,
                    s.get_impactPtr(s.impactSize,
                                    impact::impactIndices::impactVelocity));
        std::copy_n(s.get_impactPtr(0, impact::impactIndices::rawPenetration),
                    vSize - 1,
                    s.get_impactPtr(s.impactSize,
                                    impact::impactIndices::rawPenetration));
        /*std::cout<<vSize<<"\n";
        for(int i=0; i<s.impactSizeAligned; i++){
            for(int j=0; j<impact::maxColumns; j++){
                std::cout<<s.get_impact(i, j)<<" ";
            }
            std::cout<<"\n";
        }*/
        double inclination_R = M_PI / 180 * inclination;
        uint32_t length = ceil(static_cast<double>(s.postPenSize) / vSize);
        uint32_t assigned = assignThreadNum(length, nThreads);
        mtFunctionRunner(assigned, length, s.postPenSize, this,
                         &shellCalc::multiPostPen<changeDirection, fast>,
                         thickness, inclination_R, std::ref(s));

        s.completedPostPen = true;
    }
};

}  // namespace shell
