#include "shellCPP.cpp"

int main(){
    //shell test(780, .460, 2574, 1460, 6, .033, .292, 76, "Yamato");
    //shell test(780, .460, 2574, 1460, 6, .292, "Yamato");
    shell test(780, .460, 2574, 1460, 6, .292, "Yamato", 76, .033);
    //std::chrono::microseconds t1, t2;
    double total = 0;

    unsigned int runs = 1;
    for(int i=0; i<runs; i++){
        auto t1 = std::chrono::high_resolution_clock::now();
        test.calculateStd();
        auto t2 = std::chrono::high_resolution_clock::now();
        total += (double)std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
        //std::cout << duration;
    }
    test.printStdData();
    std::cout << "completed" << std::endl;
    std::cout << total / runs / 1000000000;
    //test.calculateStd();

    std::vector<double> angle;
    angle.push_back(0);
    angle.push_back(10);
    angle.push_back(20);
    angle.push_back(30);
    angle.push_back(40);
    angle.push_back(50);
    angle.push_back(60);

    /*
    for(unsigned int i=0; i<angle.size(); i++){
        printf("%f\n", angle[i]);
    }*/

    test.setAngles(&angle);
    test.calculatePostPen(400);
    test.printPostPen();
    //test.printTrajectory(2499);

}