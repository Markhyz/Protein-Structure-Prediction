#ifndef GUARD_H_EVO_ALG_MUTATION
#define GUARD_H_EVO_ALG_MUTATION

#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace mutator {
        real_chromosome_t polynomial(real_chromosome_t& chromosome, double pr, std::vector<std::pair<int, int>> bounds,
                                     uint32_t n);
        binary_chromosome_t bitFlip(binary_chromosome_t& chromosome, double pr);

    }
}

#endif