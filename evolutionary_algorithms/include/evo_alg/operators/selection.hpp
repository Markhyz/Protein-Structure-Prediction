#ifndef GUARD_H_EVO_ALG_SELECTION
#define GUARD_H_EVO_ALG_SELECTION

#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace selector {
        size_t tournament(std::vector<double> individuals_fit, size_t size = 2);
        size_t roulette(std::vector<double> individuals_fit);
    }
}

#endif