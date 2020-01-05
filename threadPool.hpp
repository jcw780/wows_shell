#include <vector>
#include <thread>
#include <functional>
#include <atomic>
#include <iostream>
#include "concurrentqueue/concurrentqueue.h" //https://github.com/cameron314/concurrentqueue

class join_threads{
    std::vector<std::thread>& threads;
    public:
    explicit join_threads(std::vector<std::thread>& threads_):
        threads(threads_) {}
    ~join_threads(){
        for(unsigned long i=0; i<threads.size(); i++){
            if(threads[i].joinable()){
                threads[i].join();
            }
        }
    }
};

template<typename T, typename A>
class parallelForEach{
    private:
    bool schedule; 
    int assigned, length;
    std::atomic<int> done;
    std::vector<std::thread> threads; 
    join_threads joiner;
    A commonArgs; //got lazy and decided to use a struct instead
    T threadFunction;
    moodycamel::ConcurrentQueue<int> work_queue;  

    void kickoffStatic(){
        for(int i=0; i<(assigned - 1); i++){
            threads.emplace_back([=] { workFunctionStatic(i); });
        } 
        workFunctionStatic(assigned - 1);
        while(done < assigned){
            std::this_thread::yield();
        }   
    }
    void workFunctionStatic(int threadID){
        //std::cout<<"start ends"<<(length*threadID/assigned)<<" "<<(length*(threadID + 1)/assigned)<<"\n";
        for(int i=(length*threadID/assigned); i<(length*(threadID + 1)/assigned); i++){
            threadFunction(i, commonArgs);
        }
        done.fetch_add(1, std::memory_order_relaxed);
    }

    void kickoffDynamic(){
        for(int i=0; i<(assigned - 1); i++){
            threads.emplace_back([=] { workFunctionDynamic(); });
        } 

        for(int i=0; i<length; i++){
            work_queue.enqueue(i);
        }
        workFunctionDynamic();
        while(done < length){
            std::this_thread::yield();
        }
        std::cout<<work_queue.size_approx()<<"\n";
    }

    void workFunctionDynamic(){
        while(done < length){
            int index;
            if(work_queue.try_dequeue(index)){
                threadFunction(index, commonArgs);
                done.fetch_add(1, std::memory_order_relaxed);
            }else{
                std::this_thread::yield();
            }
        }
    }

    public:
    parallelForEach(int threadCount, bool schedule, int length, T threadFunction, A args):joiner(threads){
        this->length = length;
        this->threadFunction = threadFunction;
        this->schedule = schedule;
        if(length > threadCount){
            assigned = threadCount;
        }else{
            assigned = length;
        }
        commonArgs = args;

        if(schedule){
            kickoffDynamic();
        }else{
            kickoffStatic();
        }
    }

    ~parallelForEach(){      
        if(schedule){
            done = length;
        }else{
            done = assigned;
        } 
    } 
};



