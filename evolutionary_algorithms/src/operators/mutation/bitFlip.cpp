#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        binary_individual_t bitFlip(binary_individual_t const& individual, double const pr) {
            binary_individual_t new_individual(individual);

            for (size_t index = 0; index < individual.getChromosome().size(); ++index) {
                double const cur_pr = utils::uniformProbGen();
                if (cur_pr < pr) {
                    new_individual.setChromosome(index, !new_individual.getChromosome()[index]);
                }
            }

            return new_individual;
        }
    }
}