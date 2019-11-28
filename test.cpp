#include "shellCPP.hpp"

#include <chrono>
#include <utility>

//Sample Test / Benchmark Function

int main(){
    //std::cout<<"Running"<<std::endl;
    //shell test(780, .460, 2574, 1460, 6, .033, .292, 76, "Yamato");
    //shell test(780, .460, 2574, 1460, 6, .292, "Yamato");
    //shell test(780, .460, 2574, 1460, 6, .292, "Yamato", 76, .033);
    //std::chrono::microseconds t1, t2;
    double total = 0.0;

    shell* test;
    shellCalc sc;
    unsigned int runs = 100;
    for(int i=0; i<runs; i++){
        test = new shell(780, .460, 2574, 1460, 6, .292, "Yamato", 76, .033);
        auto t1 = std::chrono::high_resolution_clock::now();
        //test->calculateStd();
        //std::cout<<"Running"<<std::endl;
        sc.calculateStd(*test);
        auto t2 = std::chrono::high_resolution_clock::now();
        total += (double)std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        //std::cout << duration;
        //std::cout<<"Running"<<std::endl;
        if(i < runs - 1){
            delete test;
        }
    }
    test->printStdData();
    std::cout << "completed" << std::endl;
    std::cout << total / runs / 1000000000 << std::endl;
    //test.calculateStd();

    
    std::vector<double> angle = {0, 10, 20, 30, 40, 50, 60};
    
    test->angles = std::move(angle);
    //std::cout << "completed" << std::endl;
    sc.calculatePostPen(100, *test);
    //std::cout << "completed" << std::endl;
    test->printPostPen();

    /*
    for(unsigned int i=0; i<angle.size(); i++){
        printf("%f\n", angle[i]);
    }*/

    //test.setAngles(angle);
    //test.calculatePostPen(100);
    //test.printPostPen();

    //auto var = test.postPenDataCopy();
    //std::cout<<&var<< std::endl;
    //test.printTrajectory(100);
    return 0;

}