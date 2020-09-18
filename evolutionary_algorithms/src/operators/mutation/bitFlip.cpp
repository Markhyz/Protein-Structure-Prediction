#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        binary_individual_t bitFlip(binary_individual_t const& individual, double const pr) {
            binary_individual_t mutated_individual(individual);

            binary_chromosome_t mutated_chromosome = mutated_individual.getChromosome();
            for (size_t index = 0; index < mutated_chromosome.size(); ++index) {
                double const cur_pr = utils::uniformProbGen();
                if (cur_pr < pr) {
                    mutated_chromosome[index] = !mutated_chromosome[index];
                }
            }
            mutated_individual.setChromosome(mutated_chromosome);

            return mutated_individual;
        }
    }
}