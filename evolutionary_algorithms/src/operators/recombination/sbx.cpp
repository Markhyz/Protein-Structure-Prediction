#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<real_individual_t, real_individual_t> sbx(real_individual_t& parent_1, real_individual_t& parent_2,
                                                            uint32_t n, double pr) {
            size_t const chromosome_size = parent_1.getChromosome().size();
            real_individual_t child_1(parent_1), child_2(parent_2);

            for (size_t index = 0; index < chromosome_size; ++index) {
                double const cur_pr = utils::uniform_prob_gen();
                if (cur_pr < pr) {
                    double const u = utils::uniform_prob_gen();
                    double const spread = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) : pow(1.0 / (2 * (1 - u)), 1.0 / (n + 1));
                    child_1.setChromosome(index, 0.5 * ((1 + spread) * parent_1.getChromosome()[index] +
                                                        (1 - spread) * parent_2.getChromosome()[index]));
                    child_2.setChromosome(index, 0.5 * ((1 - spread) * parent_1.getChromosome()[index] +
                                                        (1 + spread) * parent_2.getChromosome()[index]));
                }
            }

            return {child_1, child_2};
        }
    }
}