#ifndef GUARD_H_EVO_ALG_INITIALIZATION
#define GUARD_H_EVO_ALG_INITIALIZATION

#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace initializator {
        template <typename GeneType>
        using initialization_function_t = std::function<void(
            Population<Individual<GeneType>>&, typename FitnessFunction<GeneType>::const_shared_ptr, size_t const)>;

        template <typename GeneType>
        void uniformRandomInit(Population<Individual<GeneType>>& population,
                               typename FitnessFunction<GeneType>::const_shared_ptr fitness, size_t const size) {
            size_t const chromosome_size = fitness->getBounds().size();
            population.resize(size);
            for (size_t index = 0; index < size; ++index) {
                std::vector<GeneType> new_chromosome(chromosome_size);

                for (size_t gene_index = 0; gene_index < chromosome_size; ++gene_index) {
                    if constexpr (std::is_integral_v<GeneType>) {
                        std::uniform_int_distribution<GeneType> geneGenerator(
                            population.getBounds()[gene_index].first, population.getBounds()[gene_index].second);
                        new_chromosome[gene_index] = geneGenerator(utils::rng);
                    } else {
                        std::uniform_real_distribution<GeneType> geneGenerator(
                            population.getBounds()[gene_index].first, population.getBounds()[gene_index].second);
                        new_chromosome[gene_index] = geneGenerator(utils::rng);
                    }
                }

                population[index] = {fitness, new_chromosome};
            }
        }
    }
}

#endif