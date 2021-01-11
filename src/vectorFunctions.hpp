#include <immintrin.h>

#include <bitset>
#include <iostream>
#include <limits>

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdint>
#include <iomanip>

namespace wows_shell {
namespace vectorFunctions {
#ifdef __FMA__
__m128d mad(__m128d x, __m128d y, __m128d z) { return _mm_fmadd_pd(x, y, z); }
__m128d msub(__m128d x, __m128d y, __m128d z) { return _mm_fmsub_pd(x, y, z); }
#else
__m128d mad(__m128d x, __m128d y, __m128d z) {
    return _mm_add_pd(_mm_mul_pd(x, y), z);
}
__m128d msub(__m128d x, __m128d y, __m128d z) {
    return _mm_sub_pd(_mm_mul_pd(x, y), z);
}
#endif

__m128d exp2(__m128d x) {
    __m128d round_x =
        _mm_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    x = _mm_sub_pd(x, round_x);

    __m128d r = _mm_set1_pd(4.45623165388261696886670014471e-10);
    r = mad(r, x, _mm_set1_pd(7.0733589360775271430968224806e-9));
    r = mad(r, x, _mm_set1_pd(1.01780540270960163558119510246e-7));
    r = mad(r, x, _mm_set1_pd(1.3215437348041505269462510712e-6));
    r = mad(r, x, _mm_set1_pd(0.000015252733849766201174247690629));
    r = mad(r, x, _mm_set1_pd(0.000154035304541242555115696403795));
    r = mad(r, x, _mm_set1_pd(0.00133335581463968601407096905671));
    r = mad(r, x, _mm_set1_pd(0.0096181291075949686712855561931));
    r = mad(r, x, _mm_set1_pd(0.055504108664821672870565883052));
    r = mad(r, x, _mm_set1_pd(0.240226506959101382690753994082));
    r = mad(r, x, _mm_set1_pd(0.69314718055994530864272481773));
    r = mad(r, x, _mm_set1_pd(0.9999999999999999978508676375));

    constexpr auto NLD = std::numeric_limits<double>();
    const auto mantissa_bits = NLD.digits - 1;
    const double exponent_factor = (int64_t)1 << mantissa_bits;
    const auto exponent_offset_i = 2 - NLD.min_exponent;
    const double exponent_offset = (int64_t)exponent_offset_i << mantissa_bits;
    __m128d exponent = mad(round_x, _mm_set1_pd(exponent_factor),
                           _mm_set1_pd(exponent_offset));
    __m128i exponent_i = _mm_set_epi64x(static_cast<int64_t>(exponent[1]),
                                        static_cast<int64_t>(exponent[0]));

    __m128d scale = _mm_castsi128_pd(exponent_i);
    r = _mm_mul_pd(r, scale);
    return r;
}

__m128i ilogb(__m128d x) {
    constexpr int64_t exponent_mask =
        0b0111111111110000000000000000000000000000000000000000000000000000;
    __m128i exponent_raw =
        _mm_and_si128(_mm_castpd_si128(x), _mm_set1_epi64x(exponent_mask));
    __m128i exponent =
        _mm_sub_epi64(_mm_srli_epi64(exponent_raw, 52), _mm_set1_epi64x(1023));

    return exponent;
}

__m128d ldexp(const __m128d x, const __m128i exponent) {
    constexpr int64_t mantissa_mask =
        0b0000000000001111111111111111111111111111111111111111111111111111;
    __m128i exponent_o = ilogb(x);
    __m128i n_exp = _mm_sub_epi64(exponent_o, exponent);
    __m128i n_exp_r = _mm_add_epi64(n_exp, _mm_set1_epi64x(1023));

    const __m128i shifted = _mm_slli_epi64(n_exp_r, 52);

    return _mm_castsi128_pd(_mm_or_si128(
        _mm_and_si128(_mm_castpd_si128(x), _mm_set1_epi64x(mantissa_mask)),
        shifted));
}

__m128d log2(const __m128d x) {
    __m128i exponent =
        ilogb(_mm_mul_pd(x, _mm_set1_pd(1.41421356237309504880)));
    __m128d mantissa = ldexp(x, exponent);

    __m128d y = _mm_div_pd(_mm_sub_pd(mantissa, _mm_set1_pd(1)),
                           _mm_add_pd(mantissa, _mm_set1_pd(1)));
    __m128d y2 = _mm_mul_pd(y, y);
    __m128d r = _mm_set1_pd(0.243683403415639178527756320773);
    r = mad(r, y2, _mm_set1_pd(0.26136626803870009948502658));
    r = mad(r, y2, _mm_set1_pd(0.320619429891299265439389));
    r = mad(r, y2, _mm_set1_pd(0.4121983452028499242926));
    r = mad(r, y2, _mm_set1_pd(0.577078017761894161436));
    r = mad(r, y2, _mm_set1_pd(0.96179669392233355927));
    r = mad(r, y2, _mm_set1_pd(2.8853900817779295236));

    __m128d ilogbd =
        _mm_set_pd(static_cast<double>(_mm_extract_epi64(exponent, 1)),
                   static_cast<double>(_mm_extract_epi64(exponent, 0)));
    return mad(r, y, ilogbd);
}

__m128d pow(const __m128d x, const __m128d n) {
    return exp2(_mm_mul_pd(log2(x), n));
}

#ifdef __AVX2__
#ifdef __FMA__
__m256d mad(__m256d x, __m256d y, __m256d z) {
    return _mm256_fmadd_pd(x, y, z);
}
__m256d msub(__m256d x, __m256d y, __m256d z) {
    return _mm256_fmsub_pd(x, y, z);
}
#else
__m256d mad(__m256d x, __m256d y, __m256d z) {
    return _mm256_add_pd(_mm256_mul_pd(x, y), z);
}
__m256d msub(__m256d x, __m256d y, __m256d z) {
    return _mm256_sub_pd(_mm256_mul_pd(x, y), z);
}
#endif

__m256d exp2(__m256d x) {
    __m256d round_x =
        _mm256_round_pd(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    x = _mm256_sub_pd(x, round_x);

    __m256d r = _mm256_set1_pd(4.45623165388261696886670014471e-10);
    r = mad(r, x, _mm256_set1_pd(7.0733589360775271430968224806e-9));
    r = mad(r, x, _mm256_set1_pd(1.01780540270960163558119510246e-7));
    r = mad(r, x, _mm256_set1_pd(1.3215437348041505269462510712e-6));
    r = mad(r, x, _mm256_set1_pd(0.000015252733849766201174247690629));
    r = mad(r, x, _mm256_set1_pd(0.000154035304541242555115696403795));
    r = mad(r, x, _mm256_set1_pd(0.00133335581463968601407096905671));
    r = mad(r, x, _mm256_set1_pd(0.0096181291075949686712855561931));
    r = mad(r, x, _mm256_set1_pd(0.055504108664821672870565883052));
    r = mad(r, x, _mm256_set1_pd(0.240226506959101382690753994082));
    r = mad(r, x, _mm256_set1_pd(0.69314718055994530864272481773));
    r = mad(r, x, _mm256_set1_pd(0.9999999999999999978508676375));

    constexpr auto NLD = std::numeric_limits<double>();
    const auto mantissa_bits = NLD.digits - 1;
    const double exponent_factor = (int64_t)1 << mantissa_bits;
    const auto exponent_offset_i = 2 - NLD.min_exponent;
    const double exponent_offset = (int64_t)exponent_offset_i << mantissa_bits;
    __m256d exponent = mad(round_x, _mm256_set1_pd(exponent_factor),
                           _mm256_set1_pd(exponent_offset));
    __m256i exponent_i = _mm256_set_epi64x(
        static_cast<int64_t>(exponent[3]), static_cast<int64_t>(exponent[2]),
        static_cast<int64_t>(exponent[1]), static_cast<int64_t>(exponent[0]));

    __m256d scale = _mm256_castsi256_pd(exponent_i);
    r = _mm256_mul_pd(r, scale);
    return r;
}

__m256i ilogb(__m256d x) {
    constexpr int64_t exponent_mask =
        0b0111111111110000000000000000000000000000000000000000000000000000;
    __m256i exponent_raw = _mm256_and_si256(_mm256_castpd_si256(x),
                                            _mm256_set1_epi64x(exponent_mask));
    __m256i exponent = _mm256_sub_epi64(_mm256_srli_epi64(exponent_raw, 52),
                                        _mm256_set1_epi64x(1023));

    return exponent;
}

__m256d ldexp(const __m256d x, const __m256i exponent) {
    constexpr int64_t mantissa_mask =
        0b0000000000001111111111111111111111111111111111111111111111111111;
    __m256i exponent_o = ilogb(x);
    __m256i n_exp = _mm256_sub_epi64(exponent_o, exponent);
    __m256i n_exp_r = _mm256_add_epi64(n_exp, _mm256_set1_epi64x(1023));

    const __m256i shifted = _mm256_slli_epi64(n_exp_r, 52);

    return _mm256_castsi256_pd(
        _mm256_or_si256(_mm256_and_si256(_mm256_castpd_si256(x),
                                         _mm256_set1_epi64x(mantissa_mask)),
                        shifted));
}

__m256d log2(const __m256d x) {
    __m256i exponent =
        ilogb(_mm256_mul_pd(x, _mm256_set1_pd(1.41421356237309504880)));
    __m256d mantissa = ldexp(x, exponent);

    __m256d y = _mm256_div_pd(_mm256_sub_pd(mantissa, _mm256_set1_pd(1)),
                              _mm256_add_pd(mantissa, _mm256_set1_pd(1)));
    __m256d y2 = _mm256_mul_pd(y, y);
    __m256d r = _mm256_set1_pd(0.243683403415639178527756320773);
    r = mad(r, y2, _mm256_set1_pd(0.26136626803870009948502658));
    r = mad(r, y2, _mm256_set1_pd(0.320619429891299265439389));
    r = mad(r, y2, _mm256_set1_pd(0.4121983452028499242926));
    r = mad(r, y2, _mm256_set1_pd(0.577078017761894161436));
    r = mad(r, y2, _mm256_set1_pd(0.96179669392233355927));
    r = mad(r, y2, _mm256_set1_pd(2.8853900817779295236));

    __m256d ilogbd =
        _mm256_set_pd(static_cast<double>(_mm256_extract_epi64(exponent, 3)),
                      static_cast<double>(_mm256_extract_epi64(exponent, 2)),
                      static_cast<double>(_mm256_extract_epi64(exponent, 1)),
                      static_cast<double>(_mm256_extract_epi64(exponent, 0)));
    return mad(r, y, ilogbd);
}

__m256d pow(const __m256d x, const __m256d n) {
    return exp2(_mm256_mul_pd(log2(x), n));
}
#endif

}  // namespace vectorFunctions
}  // namespace wows_shell