#include <chrono>
#include <memory>
#include <utility>

#include "../shellCPP.hpp"

struct function_runtimes {
    std::chrono::nanoseconds impact;
    std::chrono::nanoseconds angle;
    std::chrono::nanoseconds post;
    std::chrono::nanoseconds dispersion;
};

template <typename Function>
auto get_runtime(Function f) {
    auto t1 = std::chrono::high_resolution_clock::now();
    f();
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
};

// Sample Test / Benchmark Function
function_runtimes runtime() {
    std::unique_ptr<wows_shell::shell> test, test1;
    wows_shell::shellCalc sc;
    sc.set_max(90.0);

    wows_shell::shellParams sp = {.460, 780, .292, 1460, 2574, 6,
                                  .033, 76,  45,   60,   0};
    wows_shell::dispersionParams dp = {10,  2.8, 1000, 5000,  0.5,
                                       0.2, 0.6, 0.8,  26630, 2.1};
    // test = std::make_unique<wows_shell::shell>(sp, dp, "Yamato");

    wows_shell::shellParams sp1 = {.457, 800, .256, 1373, 2500, 6,
                                   .033, 76,  45,   60,   0};
    wows_shell::dispersionParams dp1 = {13,   1.1, 1000, 5000,  0.6,
                                        0.25, 0.4, 0.75, 20680, 1.8};
    test = std::make_unique<wows_shell::shell>(sp1, dp1, "Kremlin");
    std::cout << wows_shell::generateHash(*test) << "\n";

    function_runtimes r;

    r.impact = get_runtime([&] {
        sc.calculateImpact<wows_shell::numerical::forwardEuler>(*test, true);
    });

    std::vector<double> angle = {0, 5, 10};

    r.post = get_runtime(
        [&] { sc.calculatePostPen(70, 0, *test, angle, true, true); });

    r.angle = get_runtime([&] { sc.calculateAngles(70, 0, *test); });

    r.dispersion = get_runtime([&] {
        sc.calculateDispersion(wows_shell::dispersion::verticalTypes::normal,
                               *test);
    });

    auto maxDist = test->maxDist();
    if (std::get<0>(maxDist) != std::numeric_limits<std::size_t>::max()) {
        std::cout << std::get<0>(maxDist) << " " << std::get<1>(maxDist)
                  << "\n";
    } else {
        std::cout << std::get<0>(maxDist) << " Error\n";
    }
    std::cout << test->interpolateDistanceImpact(
                     30000, wows_shell::impact::impactIndices::rawPenetration)
              << "\n";
    test->printTrajectory(50);
    // test->printImpactData();
    // test->printPostPenData();
    // test->printAngleData();
    // test->printDispersionData();

    return r;
}

int main() {
    long long impactTotal = 0, angleTotal = 0, postTotal = 0,
              dispersionTotal = 0;
    constexpr std::size_t runs = 1;
    for (int i = 0; i < runs; ++i) {
        auto [impact, angle, post, dispersion] = runtime();
        impactTotal += impact.count();
        angleTotal += angle.count();
        postTotal += post.count();
        dispersionTotal += dispersion.count();
        std::cout << "Stage: " << i << " Finished \n";
    }

    constexpr auto runs_d = static_cast<double>(runs);
    std::cout << "Runtimes ns\n";
    std::cout << "Impact: " << impactTotal / runs_d << "\n";
    std::cout << "Angle: " << angleTotal / runs_d << "\n";
    std::cout << "Post: " << postTotal / runs_d << "\n";
    std::cout << "Dispersion: " << dispersionTotal / runs_d << "\n";
    return 0;
}