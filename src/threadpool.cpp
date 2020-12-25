#include <atomic>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <type_traits>
#include <vector>

namespace wows_shell {
namespace utility {

class threadPool {
   private:
    std::vector<std::thread> threads;
    std::function<void(std::size_t)> f;

    std::condition_variable cv;
    std::mutex m_;
    bool ready, stop = false;

    void workerFunction(std::size_t i) {
        for (;;) {
            {
                std::unique_lock<std::mutex> lk(m_);
                cv.wait(lk, [&] { return ready || stop; });
            }
            if (stop) break;
            f(i);
        }
    }

   public:
    threadPool(std::size_t numThreads = std::thread::hardware_concurrency()) {
        threads.reserve(numThreads - 1);
        for (std::size_t i = 1; i < numThreads; ++i)
            threads.emplace_back([&, i]() { workerFunction(i); });
    }

    template <typename ThreadFunction>
    void start(ThreadFunction tf) {
        static_assert(std::is_invocable_v<ThreadFunction, std::size_t>,
                      "Function has an incorrect signature - requires "
                      "void(std::size_t).");
        f = tf;
        {
            std::lock_guard<std::mutex> lk(m_);
            ready = true;
        }
        cv.notify_all();
        tf(0);
        ready = false;
    }

    ~threadPool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop = true;
        }
        cv.notify_all();
        for (auto& t : threads) {
            t.join();
        }
    }
};

void testFunction(std::atomic<int>& counter) {
    // std::cout << "Entered \n";
    constexpr int target = 100;
    while (counter < target) {
        int v = counter.fetch_add(1, std::memory_order::memory_order_relaxed);
        if (v < target) {
            std::cout << "Run #" << v << "\n";
        }
    }
}

void testFunction2(std::atomic<int>& counter) {
    // std::cout << "Entered \n";
    constexpr int target = 200;
    while (counter < target) {
        int v = counter.fetch_add(1, std::memory_order::memory_order_relaxed);
        if (v < target) {
            std::cout << "Run #" << v << "\n";
        }
    }
}

int main() {
    std::atomic<int> counter(0);
    threadPool tp;
    tp.start([&](std::size_t i) { testFunction(counter); });

    std::cout << counter << "\n";

    tp.start([&](std::size_t i) { testFunction2(counter); });

    std::cout << counter << "\n";
    return 0;
}

}  // namespace utility
}  // namespace wows_shell