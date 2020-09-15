#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<real_individual_t, real_individual_t> arithmetic(real_individual_t const& parent_1,
                                                                   real_individual_t const& parent_2) {
            real_individual_t child_1(parent_1), child_2(parent_2);

            real_chromosome_t parent_1_chromosome = parent_1.getChromosome();
            real_chromosome_t parent_2_chromosome = parent_2.getChromosome();
            real_chromosome_t child_1_chromosome = child_1.getChromosome();
            real_chromosome_t child_2_chromosome = child_2.getChromosome();

            std::vector<std::pair<real_gene_t, real_gene_t>> bounds = parent_1.getFitnessFunction()->getBounds();

            size_t const chromosome_size = parent_1.getChromosome().size();
            double const alpha = utils::uniformProbGen();
            for (size_t index = 0; index < chromosome_size; ++index) {
                child_1_chromosome[index] =
                    parent_1_chromosome[index] * alpha + parent_2_chromosome[index] * (1 - alpha);
                child_2_chromosome[index] =
                    parent_1_chromosome[index] * (1 - alpha) + parent_2_chromosome[index] * alpha;
            }
            child_1.setChromosome(child_1_chromosome);
            child_2.setChromosome(child_2_chromosome);

            return {child_1, child_2};
        }
    }
}