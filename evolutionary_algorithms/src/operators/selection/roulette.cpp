#include "../../../include/evo_alg/commons/macros.hpp"
#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/operators/selection.hpp"

namespace evo_alg {
    namespace selector {
        size_t roulette(std::vector<double> individuals_fit) {
            std::vector<double> corrected_fit(individuals_fit);

            double const min_fit = *std::min_element(individuals_fit.begin(), individuals_fit.end());
            if (min_fit < utils::eps)
                for (double& fit : corrected_fit)
                    fit += min_fit + utils::eps;

            double const total_fit = std::accumulate(corrected_fit.begin(), corrected_fit.end(), 0.0);

            std::vector<double> normalized_fit(corrected_fit);
            for (double& fit : normalized_fit) {
                assert(fit >= utils::eps);
                fit /= total_fit;
            }

            std::vector<size_t> indexes(individuals_fit.size());
            std::iota(indexes.begin(), indexes.end(), 0);
            std::shuffle(indexes.begin(), indexes.end(), utils::rng);

            double pr = utils::uniform_prob_gen();
            double cur_pr = 0;
            size_t selected_individual = individuals_fit.size();
            for (size_t const index : indexes) {
                cur_pr += normalized_fit[index];
                if (pr - cur_pr < utils::eps) {
                    selected_individual = index;
                    break;
                }
            }
            assert(selected_individual < individuals_fit.size());

            return selected_individual;
        }
    }
}