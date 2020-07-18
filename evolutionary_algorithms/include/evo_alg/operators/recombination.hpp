#ifndef GUARD_H_EVO_ALG_RECOMBINATION
#define GUARD_H_EVO_ALG_RECOMBINATION

#include "../../../include/evo_alg/core.hpp"

#include <utility>

namespace evo_alg {
    namespace recombinator {
        template <class IndividualType>
        using crossover_function_t = std::function<std::pair<std::vector<IndividualType>, std::vector<IndividualType>>(
            std::vector<IndividualType> const&, std::vector<IndividualType> const&)>;

        std::pair<real_individual_t, real_individual_t> sbx(real_individual_t const& parent_1,
                                                            real_individual_t const& parent_2, uint32_t const n = 2,
                                                            double const pr = 0.5);
        std::pair<binary_individual_t, binary_individual_t> onePoint(binary_individual_t const& parent_1,
                                                                     binary_individual_t const& parent_2);

    }
}

#endif