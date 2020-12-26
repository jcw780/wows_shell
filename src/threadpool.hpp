#include <cassert>
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

   public:
    threadPool(std::size_t numThreads = std::thread::hardware_concurrency()) {
        threads.reserve(numThreads - 1);
        for (std::size_t i = 1; i < numThreads; ++i) {
            threads.emplace_back([&, i]() {
                for (;;) {
                    {
                        std::unique_lock<std::mutex> lk(m_);
                        cv.wait(lk, [&] { return ready || stop; });
                    }
                    if (stop) return;
                    if (ready && f) f(i);
                }
            });
        }
    }

    template <typename ThreadFunction>
    void start(ThreadFunction tf) {
        static_assert(std::is_invocable_v<ThreadFunction, std::size_t>,
                      "Function has an incorrect signature - requires "
                      "void(std::size_t).");
        {
            std::lock_guard<std::mutex> lk(m_);
            f = tf, ready = true;
        }
        cv.notify_all();
        tf(0);  // utilize main thread
        {
            std::lock_guard<std::mutex> lk(m_);
            ready = false;
        }
    }

    ~threadPool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            f = [](const std::size_t i) {};
            // pass in an empty function just in case
            stop = true, ready = false;
        }
        cv.notify_all();
        for (auto& t : threads) t.join();
    }
};
}  // namespace utility
}  // namespace wows_shell