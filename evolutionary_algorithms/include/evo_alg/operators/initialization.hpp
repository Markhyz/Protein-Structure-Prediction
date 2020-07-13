#ifndef GUARD_H_EVO_ALG_INITIALIZATION
#define GUARD_H_EVO_ALG_INITIALIZATION

#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace initializator {
        template <typename GeneType>
        std::vector<std::vector<GeneType>> uniformRandomInit(std::vector<std::vector<GeneType>> chromosomes,
                                                             std::vector<std::pair<GeneType, GeneType>> bounds,
                                                             size_t size) {
            chromosomes.resize(size);
            if constexpr (std::is_integral_v<GeneType>) {
                for (size_t index = 0; index < size; ++index) {
                    std::uniform_int_distribution<GeneType> geneGenerator(bounds[index].first, bounds[index].second);
                    chromosomes[i] = geneGenerator(utils::rng);
                }
            } else {
                for (size_t index = 0; index < size; ++index) {
                    std::uniform_real_distribution<GeneType> geneGenerator(bounds[index].first, bounds[index].second);
                    chromosomes[i] = geneGenerator(utils::rng);
                }
            }
        }
    }
}

#endif