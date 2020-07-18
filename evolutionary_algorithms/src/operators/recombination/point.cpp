#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<binary_individual_t, binary_individual_t> onePoint(binary_individual_t const& parent_1,
                                                                     binary_individual_t const& parent_2) {
            size_t const chromosome_size = parent_1.getChromosome().size();
            binary_individual_t child_1(parent_1), child_2(parent_2);
            std::uniform_int_distribution<uint32_t> point_dist(1, (uint32_t) chromosome_size - 1);

            uint32_t const cut_point = point_dist(utils::rng);
            for (size_t index = cut_point; index < chromosome_size; ++index) {
                child_1.setChromosome(index, parent_2.getChromosome()[index]);
                child_2.setChromosome(index, parent_1.getChromosome()[index]);
            }

            return {child_1, child_2};
        }
    }
}