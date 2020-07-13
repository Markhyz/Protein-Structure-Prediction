#include "../../../include/evo_alg/commons/macros.hpp"
#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/mutation.hpp"

namespace evo_alg {
    namespace mutator {
        real_chromosome_t polynomial(real_chromosome_t& chromosome, double pr, std::vector<std::pair<int, int>> bounds,
                                     uint32_t n) {
            real_chromosome_t new_chromosome(chromosome);

            for (size_t index = 0; index < chromosome.size(); ++index) {
                double const cur_pr = utils::uniform_prob_gen();
                if (cur_pr < pr) {
                    double const u = utils::uniform_prob_gen();
                    double const delta = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) - 1 : 1 - pow((2 * (1 - u)), 1 / (n + 1));
                    new_chromosome[index] += delta * (u < 0.5 ? chromosome[index] - bounds[index].first
                                                              : bounds[index].second - chromosome[index]);
                }
                assert(new_chromosome[index] - bounds[index].first > 0);
                assert(bounds[index].second - new_chromosome[index] > 0);
            }

            return new_chromosome;
        }
    }
}