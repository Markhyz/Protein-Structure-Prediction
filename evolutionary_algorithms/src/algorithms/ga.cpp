#include "../../include/evo_alg/commons/utils.hpp"
#include "../../include/evo_alg/core.hpp"
#include "../../include/evo_alg/operators.hpp"

namespace evo_alg {
    template <typename... ChromosomeTypes>
    Population<Individual<ChromosomeTypes...>>
    ga(size_t const iterations, size_t const pop_size,
       typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr const fitness,
       typename recombinator::crossover_function_t<Individual<ChromosomeTypes...>> const& crossover_fun,
       double const crossover_pr,
       typename mutator::mutation_function_t<Individual<ChromosomeTypes...>> const& mutation_fun,
       double const mutation_pr, size_t const elite_size,
       typename initializator::initialization_function_t<Individual<ChromosomeTypes...>> const& init_fun,
       typename selector::selection_function_t<Individual<ChromosomeTypes...>> const& selection_fun) {
        Population<Individual<ChromosomeTypes...>> population;

        init_fun(population, fitness, pop_size);
        population.evaluateFitness();

        for (size_t it = 0; it < iterations; ++it) {
            std::vector<double> pop_fitness;
            for (size_t index = 0; index < pop_size; ++index)
                pop_fitness.push_back(population.getIndividual(index).getFitnessValue());

            Population<Individual<ChromosomeTypes...>> new_population;
            while (new_population.size() < pop_size) {
                size_t const parent_1 = selection_fun(pop_fitness);
                pop_fitness[parent_1] = -1;
                size_t const parent_2 = selection_fun(pop_fitness);

                double const pr = utils::uniformProbGen();
                Individual<ChromosomeTypes...> child_1;
                Individual<ChromosomeTypes...> child_2;
                if (pr < crossover_pr) {
                    tie(child_1, child_2) =
                        crossover_fun(population.getIndividual(parent_1), population.getIndividual(parent_2));
                } else {
                    child_1 = population.getIndividual(parent_1);
                    child_2 = population.getIndividual(parent_2);
                }

                child_1 = mutation_fun(child_1, mutation_pr);
                child_2 = mutation_fun(child_2, mutation_pr);

                new_population.appendIndividual(child_1);
                if (new_population.size() < pop_size)
                    new_population.appendIndividual(child_2);
            }

            new_population.evaluateFitness();

            if (elite_size > 0) {
                const std::vector<size_t> sorted_individuals = population.getSortedIndividuals();
                const std::vector<size_t> new_sorted_individuals = new_population.getSortedIndividuals();
                for (size_t index = 0; index < elite_size; ++index)
                    new_population[new_sorted_individuals[pop_size - index - 1]] =
                        population[sorted_individuals[index]];
            }

            population = std::move(new_population);
        }

        return population;
    }
}