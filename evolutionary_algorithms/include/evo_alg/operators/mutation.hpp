#ifndef GUARD_H_EVO_ALG_MUTATION
#define GUARD_H_EVO_ALG_MUTATION

#include "../../../include/evo_alg/core.hpp"

namespace evo_alg {
    namespace mutator {
        template <class IndividualType>
        using mutation_function_t =
            std::function<std::vector<IndividualType>(std::vector<IndividualType> const&, double const)>;

        real_individual_t polynomial(real_individual_t const& individual, double const pr, uint32_t const n = 20);
        binary_individual_t bitFlip(binary_individual_t const& individual, double const pr);

    }
}

#endif