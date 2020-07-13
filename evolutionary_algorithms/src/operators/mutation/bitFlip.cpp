#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        binary_chromosome_t bitFlip(binary_chromosome_t& chromosome, double pr) {
            binary_chromosome_t new_chromosome(chromosome);

            for (size_t index = 0; index < chromosome.size(); ++index) {
                double const cur_pr = utils::uniform_prob_gen();
                if (cur_pr < pr) {
                    new_chromosome[index] = !new_chromosome[index];
                }
            }

            return new_chromosome;
        }
    }
}