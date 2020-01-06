#include <string>
#include <cstdlib>
#include <ctime>
#include "threadPool.hpp"

typedef struct {
    int str;
    std::vector<int>* target;
}Args;


void test(int i, Args out){
    //std::cout<<i<<" "<<out.str<<"\n";
    (out.target->data())[i] = 100;
}

bool testFunc(){
    std::cout<<"entered\n";
    int testSize = rand();
    std::function<void(int, Args)> fTest = test;
    std::vector<int> testCompletion(testSize);

    Args a;
    a.str = 10;
    a.target = &testCompletion;
    std::cout<<"entered\n";
    parallelForEach<std::function<void(int, Args)>, Args> stp (12, true, testSize, test, a );
    std::cout<<"entered\n";
    bool correct = true;
    for(int i=0; i<testCompletion.size(); i++){
        if(testCompletion[i] != 100){
            correct = false;
        }
    }
    std::cout<<testSize<<" "<<correct<<"\n";
    return correct;
}

int main(){
    //std::string str = "Print";
    int total = 0;
    srand (time(NULL));
    for(int runs=0; runs<10; runs++){
        total += testFunc();
    }
    std::cout<<total<<std::endl;
    
    return 0;
}