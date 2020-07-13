#ifndef GUARD_H_EVO_ALG_RECOMBINATION
#define GUARD_H_EVO_ALG_RECOMBINATION

#include "../../../include/evo_alg/core.hpp"

#include <utility>

namespace evo_alg {
    namespace recombinator {
        std::pair<real_chromosome_t, real_chromosome_t> sbx(real_chromosome_t& parent_1, real_chromosome_t& parent_2,
                                                            uint32_t n, double pr = 0.5);
        std::pair<binary_chromosome_t, binary_chromosome_t> onePoint(binary_chromosome_t& parent_1,
                                                                     binary_chromosome_t& parent_2);

    }
}

#endif