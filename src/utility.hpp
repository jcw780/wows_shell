#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iterator>
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

std::string base64_encode(const std::vector<char> &in) {
    const std::string b = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;

    int val=0, valb=-6;
    for (char c : in) {
        val = (val<<8) + c;
        valb += 8;
        while (valb>=0) {
            out.push_back(b[(val>>valb)&0x3F]);
            valb-=6;
        }
    }
    if (valb>-6) out.push_back(b[((val<<8)>>(valb+8))&0x3F]);
    while (out.size()%4) out.push_back('=');
    return out;
}
}  // namespace utility
}  // namespace wows_shell