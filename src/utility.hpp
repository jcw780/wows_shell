#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iterator>
namespace shell {
namespace utility {
template <typename>
struct falsy {
    constexpr static auto value = false;
};
template <typename T>
constexpr inline auto falsy_v = falsy<T>::value;
}  // namespace utility
}  // namespace shell