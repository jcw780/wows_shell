#include <cmath>

namespace shell{
    namespace utility{
        template <std::size_t N>
        void fmaArrInplace(double x, double *y, double *z) {
            for (unsigned int i = 0; i < N; i++) {
                z[i] = std::fma(x, y[i], z[i]);
            }
        }

        template <std::size_t N>
        void fmaArr(double x, double *y, double *z, double *out) {
            for (unsigned int i = 0; i < N; i++) {
                out[i] = std::fma(x, y[i], z[i]);
            }
        }

        template <std::size_t N>
        void multiplyArrInplace(double x, double *z) {
            for (unsigned int i = 0; i < N; i++) {
                z[i] *= x;
            }
        }

        template <std::size_t N>
        void multiplyArr(double x, double *z, double *output) {
            for (unsigned int i = 0; i < N; i++) {
                output[i] = x * z[i];
            }
        }

        template <std::size_t N>
        void addArrInplace(double *x, double *z) {
            for (unsigned int i = 0; i < N; i++) {
                z[i] += x[i];
            }
        }

        template <std::size_t N>
        void addArr(double *x, double *z, double *output) {
            for (unsigned int i = 0; i < N; i++) {
                output[i] = x[i] + z[i];
            }
        }
    }
}