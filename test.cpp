#include "shellCPP.hpp"

#include <chrono>
#include <utility>

// Sample Test / Benchmark Function
void runtime() {
    double total = 0.0;

    shell::shell *test;
    shell::shellCalc sc;
    unsigned int runs = 10;
    for (int i = 0; i < runs; i++) {
        test = new shell::shell(.460, 780, .292, 1460, 2574, 6, .033, 76, 45,
                                60, "Yamato");
        auto t1 = std::chrono::high_resolution_clock::now();
        sc.calculateImpact<shell::numerical::rungeKutta4>(*test, true);
        auto t2 = std::chrono::high_resolution_clock::now();
        total += (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
                     t2 - t1)
                     .count();
        if (i < runs - 1) {
            delete test;
        }
    }

    // std::cout << "completed" << std::endl;
    // test.calculateStd();

    std::vector<double> angle = {0, 10};
    // std::cout<<"Started\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    sc.calculatePostPen(10, -20, *test, angle, true, true);
    auto t2 = std::chrono::high_resolution_clock::now();

    sc.calculateAngles(410, -20, *test);

    test->printImpactData();
    // test->printPostPenData();
    // test->printAngleData();

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