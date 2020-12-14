#pragma once

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <string>
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

    auto addToOutput = [&](const char* target, uint8_t N) {
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
    return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
}

double cdf(const double& x) { return (1 + erf(x / sqrt(2))) / 2; }

// https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
// This is only the single precision version but it's good enough
double MBG_erfinv(double x) {
    double p;
    double w = -log((1.0f - x) * (1.0f + x));
    if (w < 5.000000f) {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    } else {
        w = sqrtf(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

double invCDF(double x) { return sqrt(2) * MBG_erfinv(2 * x - 1); }

}  // namespace utility
}  // namespace wows_shell