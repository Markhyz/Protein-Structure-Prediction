#include "../../../include/evo_alg/commons/macros.hpp"
#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        real_individual_t polynomial(real_individual_t& individual, double pr, std::vector<std::pair<int, int>> bounds,
                                     uint32_t n) {
            real_individual_t new_individual(individual);

            for (size_t index = 0; index < individual.getChromosome().size(); ++index) {
                double const cur_pr = utils::uniform_prob_gen();
                if (cur_pr < pr) {
                    double const u = utils::uniform_prob_gen();
                    double const delta = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) - 1 : 1 - pow((2 * (1 - u)), 1 / (n + 1));
                    new_individual.setChromosome(
                        index, new_individual.getChromosome()[index] +
                                   delta * (u < 0.5 ? individual.getChromosome()[index] - bounds[index].first
                                                    : bounds[index].second - individual.getChromosome()[index]));
                }
                assert(new_individual.getChromosome()[index] - bounds[index].first > 0);
                assert(bounds[index].second - new_individual.getChromosome()[index] > 0);
            }

            return new_individual;
        }
    }
}