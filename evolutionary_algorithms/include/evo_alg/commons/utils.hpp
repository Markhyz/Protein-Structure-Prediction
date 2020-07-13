#ifndef GUARD_H_EVO_ALG_UTILS
#define GUARD_H_EVO_ALG_UTILS

#include <chrono>
#include <random>

namespace evo_alg {
    namespace utils {
        constexpr double eps = 1e-9;

        int64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        std::mt19937_64 rng(ss);
        std::uniform_real_distribution<double> uniform_dist(0, 1);
        auto uniform_prob_gen = []() { return uniform_dist(rng); };
    }
}

#endif