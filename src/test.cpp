#include "shellCPP.hpp"

#include <memory>
#include <chrono>
#include <utility>

// Sample Test / Benchmark Function
void runtime() {
    double total = 0.0;

    //shell::shell *test;
    std::unique_ptr<shell::shell> test;
    shell::shellCalc sc;
    unsigned int runs = 1;
    for (int i = 0; i < runs; i++) {

        //test = new shell::shell(.460, 780, .292, 1460, 2574, 6, .001, 2, 91, 60,
        //                        76, "Yamato");
        test = std::make_unique<shell::shell>(
            .460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0, "Yamato");


        // test = new shell::shell(.102, 805, .3536, 15.2, 2300, 10, .01, 17,
        // 45,
        //                        60, "Yamato");
        auto t1 = std::chrono::high_resolution_clock::now();
        sc.calculateImpact<shell::numerical::forwardEuler, false>(*test, true, 1);
        auto t2 = std::chrono::high_resolution_clock::now();
        total += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
                     t2 - t1)
                     .count();
    }

    // std::cout << "completed" << std::endl;
    // test.calculateStd();

    std::vector<double> angle = {0, 5, 10};
    // std::cout<<"Started\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    sc.calculatePostPen(70, 0, *test, angle, true, true);
    auto t2 = std::chrono::high_resolution_clock::now();

    sc.calculateAngles(76, 0, *test);

    //test->printImpactData();
    // std::cout << test->interpolateDistanceImpact(
    //                 30000, shell::impact::impactDataIndex::rawPen)
    //          << "\n";
    // test->printTrajectory(0);
    test->printPostPenData();
    //test->printAngleData();

    std::cout << std::fixed << std::setprecision(10)
              << total / runs / 1000000000 << std::endl;
    std::cout << std::fixed << std::setprecision(10)
              << (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
                     t2 - t1)
                         .count() /
                     1000000000
              << "\n";
    // std::cout<<"Ended\n";
}

int main() {
    runtime();
    return 0;
}