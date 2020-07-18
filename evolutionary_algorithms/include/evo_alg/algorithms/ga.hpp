#ifndef GUARD_H_EVO_ALG_GA
#define GUARD_H_EVO_ALG_GA

#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/core.hpp"
#include "../../../include/evo_alg/operators.hpp"

namespace evo_alg {
    template <typename... ChromosomeTypes>
    Population<Individual<ChromosomeTypes...>>
    ga(size_t const iterations, size_t const pop_size,
       typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr const fitness,
       typename recombinator::crossover_function_t<Individual<ChromosomeTypes...>> const& crossover_fun,
       double const crossover_pr,
       typename mutator::mutation_function_t<Individual<ChromosomeTypes...>> const& mutation_fun,
       double const mutation_pr, size_t const elite_size = 1,
       typename initializator::initialization_function_t<Individual<ChromosomeTypes...>> const& init_fun =
           initializator::uniformRandomInit,
       typename selector::selection_function_t<Individual<ChromosomeTypes...>> const& selection_fun =
           selector::tournament);
}

#endif