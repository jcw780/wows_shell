#pragma once

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
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

uint32_t reverse_bytes(uint32_t bytes)
{
    uint32_t aux = 0;
    uint8_t byte;
    int i;

    for(i = 0; i < 32; i+=8)
    {
        byte = (bytes >> i) & 0xff;
        aux |= byte << (32 - 8 - i);
    }
    return aux;
}

template <typename T>
std::string base85Encode(const T& in){
    static_assert(std::is_same<typename T::value_type, char>(), "Only accepts char elements");
    const std::string charSet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&()*+-;<=>?@^_`{|}~";
    std::string output;

    constexpr std::uint32_t Pow85[5] = {
		52200625ul, 614125ul, 7225ul, 85ul, 1ul
	};

    auto mapToCharSet = [&](uint32_t InTuple, uint8_t stage){
        return charSet[((InTuple / Pow85[stage]) % 85ul)];
    };

    auto addToOutput = [&](const char* target){
        const std::uint32_t InTuple = 
        reverse_bytes(
			*reinterpret_cast<const std::uint32_t*>(target)
		);

        for (uint8_t j=0; j<5; ++j){
            output += mapToCharSet(InTuple, j);
        }
    };

    for( uint8_t i = 0; i < in.size() / 4; ++i )
        addToOutput(&in[i * 4]);

    std::array<char, 4> temp{0};
    for( std::size_t i = in.size() - in.size() % 4; i < in.size(); ++i){
        temp[4-i] = in[i];
    }

    addToOutput(temp.data());

    return output;
}


}  // namespace utility
}  // namespace wows_shell