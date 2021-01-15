#pragma once

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

namespace wows_shell {
namespace utility {
template <typename>
struct falsy {
    constexpr static auto value = false;
};
template <typename T>
constexpr inline auto falsy_v = falsy<T>::value;

template <typename T>
std::string base64_encode(const T& in) {
    static_assert(std::is_same<typename T::value_type, char>(),
                  "Only accepts char elements");
    const std::string b =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;

    int val = 0, valb = -6;
    for (char c : in) {
        val = (val << 8) + c;
        valb += 8;
        while (valb >= 0) {
            out.push_back(b[(val >> valb) & 0x3F]);
            valb -= 6;
        }
    }
    if (valb > -6) out.push_back(b[((val << 8) >> (valb + 8)) & 0x3F]);
    while (out.size() % 4) out.push_back('=');
    return out;
}

template <typename T>
std::string base85Encode(const T& in, const bool pad = false) {
    static_assert(std::is_same<typename T::value_type, char>(),
                  "Only accepts char elements");
    const std::string charSet =
        "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&()*"
        "+-;<=>?@^_`{|}~";
    std::string output;

    constexpr uint32_t Pow85[5] = {52200625ul, 614125ul, 7225ul, 85ul, 1ul};

    auto mapToCharSet = [&](uint32_t InTuple, uint8_t stage) {
        return charSet[((InTuple / Pow85[stage]) % 85ul)];
    };

    auto addToOutput = [&](const char* target, const uint8_t N) {
        uint32_t inTuple = 0;

        for (uint8_t j = 0; j < (sizeof(uint32_t) / sizeof(char)); ++j) {
            inTuple <<= 8;
            inTuple += target[j];
        }

        for (uint8_t j = 0; j < (N + 1); ++j) {
            output += mapToCharSet(inTuple, j);
        }
    };

    std::div_t sizeBlockRemainder = std::div(in.size(), 4);
    uint8_t remainder = sizeBlockRemainder.rem;
    std::size_t blocks = sizeBlockRemainder.quot;

    for (uint8_t i = 0; i < blocks; ++i) addToOutput(&in[i * 4], 4);

    std::array<char, 4> temp{0};
    for (std::size_t i = in.size() - remainder; i < in.size(); ++i) {
        temp[4 - i] = in[i];
    }

    addToOutput(temp.data(), pad ? 4 : remainder);

    return output;
}

double pdf(const double& x) {
    return (1.0l / sqrt(2.0l * M_PI)) * exp(-0.5l * x * x);
}

double cdf(const double& x) { return (1 + erf(x / sqrt(2))) / 2; }

// https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
// This is only the single precision version but it's good enough
double MBG_erfinv(double x) {
    double p;
    double w = -log((1.0l - x) * (1.0l + x));
    if (w < 5.000000l) {
        w = w - 2.500000l;
        p = 2.81022636e-08l;
        p = 3.43273939e-07l + p * w;
        p = -3.5233877e-06l + p * w;
        p = -4.39150654e-06l + p * w;
        p = 0.00021858087l + p * w;
        p = -0.00125372503l + p * w;
        p = -0.00417768164l + p * w;
        p = 0.246640727l + p * w;
        p = 1.50140941l + p * w;
    } else {
        w = sqrt(w) - 3.000000l;
        p = -0.000200214257l;
        p = 0.000100950558l + p * w;
        p = 0.00134934322l + p * w;
        p = -0.00367342844l + p * w;
        p = 0.00573950773l + p * w;
        p = -0.0076224613l + p * w;
        p = 0.00943887047l + p * w;
        p = 1.00167406l + p * w;
        p = 2.83297682l + p * w;
    }
    return p * x;
}

double invCDF(double x) { return sqrt(2) * MBG_erfinv(2 * x - 1); }

class threadPool {
   private:
    std::atomic<std::size_t> finished = {0};
    std::vector<std::thread> threads;
    std::function<void(std::size_t)> f;

    std::condition_variable cv, cv_finished;
    std::mutex m_;
    std::atomic_size_t busy;
    bool ready = false, stop = false;

   public:
    threadPool(std::size_t numThreads = std::thread::hardware_concurrency()) {
        finished.store(0, std::memory_order_release);
        threads.reserve(numThreads - 1);
        for (std::size_t i = 1; i < numThreads; ++i) {
            threads.emplace_back([&, i]() {
                for (;;) {
                    std::unique_lock<std::mutex> lk(m_);
                    cv.wait(lk, [&] { return ready || stop; });

                    if (stop) return;
                    if (ready && f) {
                        busy.fetch_add(1, std::memory_order_relaxed);
                        lk.unlock();
                        f(i);
                        lk.lock();
                        busy.fetch_sub(1, std::memory_order_relaxed);
                        cv_finished.notify_one();
                    }
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
            finished.store(0, std::memory_order_acq_rel);
        }
        cv.notify_all();
        tf(0);  // utilize main thread
        auto fetched = finished.fetch_add(1, std::memory_order_acq_rel);
        if (fetched >= threads.size()) {
            cv.notify_all();
            finished.store(0, std::memory_order_acquire);
        } else {
            std::unique_lock<std::mutex> lk(m_);
            cv.wait(lk, [&] {
                return finished.load(std::memory_order_acquire) != fetched;
            });
        }
        std::cout << "Done\n";
        {
            std::unique_lock<std::mutex> lk(m_);
            cv_finished.wait(lk, [&] {
                return (busy.load(std::memory_order_relaxed) == 0);
            });
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