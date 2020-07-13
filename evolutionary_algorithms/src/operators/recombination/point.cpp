#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<binary_chromosome_t, binary_chromosome_t> onePoint(binary_chromosome_t& parent_1,
                                                                     binary_chromosome_t& parent_2) {
            size_t const chromosome_size = parent_1.size();
            binary_chromosome_t child_1(parent_1), child_2(parent_2);
            std::uniform_int_distribution<uint32_t> point_dist(1, (uint32_t) chromosome_size - 1);

            uint32_t const cut_point = point_dist(utils::rng);
            for (size_t index = cut_point; index < chromosome_size; ++index) {
                child_1[index] = parent_2[index];
                child_2[index] = parent_1[index];
            }

            return {child_1, child_2};
        }
    }
}