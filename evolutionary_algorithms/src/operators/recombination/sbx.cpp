#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/recombination.hpp"

namespace evo_alg {
    namespace recombinator {
        std::pair<real_individual_t, real_individual_t>
        sbx(real_individual_t const& parent_1, real_individual_t const& parent_2, uint32_t const n, double const pr) {
            real_individual_t child_1(parent_1), child_2(parent_2);

            real_chromosome_t parent_1_chromosome = parent_1.getChromosome();
            real_chromosome_t parent_2_chromosome = parent_2.getChromosome();
            real_chromosome_t child_1_chromosome = child_1.getChromosome();
            real_chromosome_t child_2_chromosome = child_2.getChromosome();

            std::vector<std::pair<real_gene_t, real_gene_t>> bounds = parent_1.getFitnessFunction()->getBounds();

            size_t const chromosome_size = parent_1.getChromosome().size();
            for (size_t index = 0; index < chromosome_size; ++index) {
                double const cur_pr = utils::uniformProbGen();
                if (cur_pr < pr) {
                    double const u = utils::uniformProbGen();
                    double const spread = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) : pow(1.0 / (2 * (1 - u)), 1.0 / (n + 1));

                    child_1_chromosome[index] =
                        0.5 * ((1 + spread) * parent_1_chromosome[index] + (1 - spread) * parent_2_chromosome[index]);
                    if (utils::numericLower(child_1_chromosome[index], bounds[index].first))
                        child_1_chromosome[index] = bounds[index].first;
                    if (utils::numericGreater(child_1_chromosome[index], bounds[index].second))
                        child_1_chromosome[index] = bounds[index].second;

                    child_2_chromosome[index] =
                        0.5 * ((1 - spread) * parent_1_chromosome[index] + (1 + spread) * parent_2_chromosome[index]);
                    if (utils::numericLower(child_2_chromosome[index], bounds[index].first))
                        child_2_chromosome[index] = bounds[index].first;
                    if (utils::numericGreater(child_2_chromosome[index], bounds[index].second))
                        child_2_chromosome[index] = bounds[index].second;
                }
            }
            child_1.setChromosome(child_1_chromosome);
            child_2.setChromosome(child_2_chromosome);

            return {child_1, child_2};
        }
    }
}