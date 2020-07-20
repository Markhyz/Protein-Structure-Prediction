#include "../../../include/evo_alg/commons/macros.hpp"
#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        real_individual_t polynomial(real_individual_t const& individual, double const pr, uint32_t const n) {
            real_individual_t mutated_individual(individual);

            real_chromosome_t individual_chromosome = individual.getChromosome();
            real_chromosome_t mutated_chromosome = mutated_individual.getChromosome();
            std::vector<std::pair<double, double>> bounds = individual.getBounds();
            for (size_t index = 0; index < mutated_chromosome.size(); ++index) {
                double const cur_pr = utils::uniformProbGen();
                if (cur_pr < pr) {
                    double const u = utils::uniformProbGen();
                    double const delta = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) - 1 : 1 - pow((2 * (1 - u)), 1 / (n + 1));
                    mutated_chromosome[index] +=
                        delta * (u < 0.5 ? individual_chromosome[index] - bounds[index].first
                                         : bounds[index].second - individual_chromosome[index]);
                }
                assert(utils::numericGreaterEqual(mutated_chromosome[index], bounds[index].first));
                assert(utils::numericLowerEqual(mutated_chromosome[index], bounds[index].second));
            }
            mutated_individual.setChromosome(mutated_chromosome);

            return mutated_individual;
        }
    }
}