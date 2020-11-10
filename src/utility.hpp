#pragma once

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <iterator>
#include <string>
#include <vector>
#include <array>

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
    static_assert(std::is_same<typename T::value_type, char>(), "Only accepts char elements");
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
  
template <typename T>
std::string base85Encode(const T& in, const bool pad=false){
    static_assert(std::is_same<typename T::value_type, char>(), "Only accepts char elements");
    const std::string charSet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&()*+-;<=>?@^_`{|}~";
    std::string output;

    constexpr uint32_t Pow85[5] = {
		52200625ul, 614125ul, 7225ul, 85ul, 1ul
	};

    auto mapToCharSet = [&](uint32_t InTuple, uint8_t stage){
        return charSet[((InTuple / Pow85[stage]) % 85ul)];
    };

    auto addToOutput = [&](const char* target, uint8_t N){
        uint32_t inTuple = 0;

        for (uint8_t j=0; j<(sizeof(uint32_t) / sizeof(char)); ++j){
            inTuple <<= 8;
            inTuple += target[j];
        }

        for (uint8_t j=0; j<(N+1); ++j){
            output += mapToCharSet(inTuple, j);
        }
    };

    std::div_t sizeBlockRemainder = std::div(in.size(), 4);
    uint8_t remainder = sizeBlockRemainder.rem;
    std::size_t blocks = sizeBlockRemainder.quot;

    for( uint8_t i = 0; i < blocks; ++i )
        addToOutput(&in[i * 4], 4);

    std::array<char, 4> temp{0};
    for( std::size_t i = in.size() - remainder; i < in.size(); ++i){
        temp[4-i] = in[i];
    }

    addToOutput(temp.data(), pad ? 4 : remainder);

    return output;
}


}  // namespace utility
}  // namespace wows_shell