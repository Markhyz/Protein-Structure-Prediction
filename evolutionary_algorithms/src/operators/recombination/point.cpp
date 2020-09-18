#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<binary_individual_t, binary_individual_t> onePoint(binary_individual_t const& parent_1,
                                                                     binary_individual_t const& parent_2) {
            binary_individual_t child_1(parent_1), child_2(parent_2);

            binary_chromosome_t parent_1_chromosome = parent_1.getChromosome();
            binary_chromosome_t parent_2_chromosome = parent_2.getChromosome();
            binary_chromosome_t child_1_chromosome = child_1.getChromosome();
            binary_chromosome_t child_2_chromosome = child_2.getChromosome();

            size_t const chromosome_size = parent_1.getChromosome().size();
            std::uniform_int_distribution<size_t> point_dist(1, chromosome_size - 1);
            size_t const cut_point = point_dist(utils::rng);
            for (size_t index = cut_point; index < chromosome_size; ++index) {
                child_1_chromosome[index] = parent_2_chromosome[index];
                child_2_chromosome[index] = parent_1_chromosome[index];
            }
            child_1.setChromosome(child_1_chromosome);
            child_2.setChromosome(child_2_chromosome);

            return {child_1, child_2};
        }
    }
}