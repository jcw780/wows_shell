#include <chrono>
#include <memory>
#include <utility>

#include "../shellCPP.hpp"

// Sample Test / Benchmark Function
void runtime() {
    double total = 0.0;

    // wows_shell::shell *test;
    std::unique_ptr<wows_shell::shell> test;
    wows_shell::shellCalc sc;
    sc.set_max(25.0);
    unsigned int runs = 1;
    test = std::make_unique<wows_shell::shell>(.460, 780, .292, 1460, 2574, 6, .033,
                                          76, 45, 60, 0, "Yamato");
    
    std::cout << wows_shell::generateHash(*test) << "\n";
    for (unsigned int i = 0; i < runs; i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        sc.calculateImpact<wows_shell::numerical::adamsBashforth5, false>(*test,
                                                                     false);
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

    test->printImpactData();
    std::size_t maxIndex = test->maxDist();
    if (maxIndex != std::numeric_limits<std::size_t>::max()) {
        std::cout << maxIndex << " "
                  << test->get_impact(maxIndex,
                                      wows_shell::impact::impactIndices::distance)
                  << "\n";
    } else {
        std::cout << maxIndex << " Error\n";
    }
    std::cout << test->interpolateDistanceImpact(
                     30000, wows_shell::impact::impactIndices::rawPenetration)
              << "\n";
    // test->printTrajectory(0);
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

    std::cout<<wows_shell::utility::base85Encode(std::string("Hello"))<<"\n";
}

int main() {
    runtime();
    return 0;
}