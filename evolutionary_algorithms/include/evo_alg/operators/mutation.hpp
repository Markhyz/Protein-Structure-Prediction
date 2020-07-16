#ifndef GUARD_H_EVO_ALG_MUTATION
#define GUARD_H_EVO_ALG_MUTATION

#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace mutator {
        real_individual_t polynomial(real_individual_t& individual, double pr, std::vector<std::pair<int, int>> bounds,
                                     uint32_t n);
        binary_individual_t bitFlip(binary_individual_t& individual, double pr);

    }
}

#endif