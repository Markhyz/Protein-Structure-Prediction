#ifndef GUARD_H_EVO_ALG_UTILS
#define GUARD_H_EVO_ALG_UTILS

#include <chrono>
#include <random>

namespace evo_alg {
    namespace utils {
        constexpr double eps = 1e-15;

        extern std::mt19937_64 rng;

        extern int mpi_size, mpi_rank;

        double uniformProbGen();

        void initMPI(int const size, int const rank);

        template <typename T>
        bool numericLower(T const x, T const y, double const precision = eps) {
            return x - y < (T) -precision;
        }

        template <typename T>
        bool numericGreater(T const x, T const y, double const precision = eps) {
            return x - y > (T) precision;
        }

        template <typename T>
        bool numericEqual(T const x, T const y, double const precision = eps) {
            return x - y <= (T) precision && x - y >= (T) -precision;
        }

        template <typename T>
        bool numericLowerEqual(T const x, T const y, double const precision = eps) {
            return numericLower(x, y, precision) || numericEqual(x, y, precision);
        }

        template <typename T>
        bool numericGreaterEqual(T const x, T const y, double const precision = eps) {
            return numericGreater(x, y, precision) || numericEqual(x, y, precision);
        }
    }
}

#endif