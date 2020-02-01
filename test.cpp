#include "shellCPP.hpp"

#include <chrono>
#include <utility>

//Sample Test / Benchmark Function

int main(){
    double total = 0.0;

    shell::shell* test;
    shell::shellCalc sc;
    unsigned int runs = 10;
    for(int i=0; i<runs; i++){
        test = new shell::shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76, .033);
        auto t1 = std::chrono::high_resolution_clock::now();
        sc.calculateImpact(*test, true);
        auto t2 = std::chrono::high_resolution_clock::now();
        total += (double)std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        if(i < runs - 1){
            delete test;
        }
    }

    //std::cout << "completed" << std::endl;
    //test.calculateStd();

    
    std::vector<double> angle = {10};
    //std::cout<<"Started\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    sc.calculatePostPen(70, -20, *test, angle, true, true);
    auto t2 = std::chrono::high_resolution_clock::now();

    test->printImpactData();
    test->printPostPenData();

    std::cout <<std::fixed<<std::setprecision(10)<< total / runs / 1000000000 << std::endl;
    std::cout<<std::fixed<<std::setprecision(10)<< (double)std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count() / 1000000000 <<"\n";
    //std::cout<<"Ended\n";

    return 0;

}