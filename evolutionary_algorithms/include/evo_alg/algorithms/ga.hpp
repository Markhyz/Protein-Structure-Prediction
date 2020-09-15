#ifndef GUARD_H_EVO_ALG_GA
#define GUARD_H_EVO_ALG_GA

#include "../../../include/evo_alg/commons/utils.hpp"
#include "../../../include/evo_alg/core.hpp"
#include "../../../include/evo_alg/operators.hpp"

#include <omp.h>

#include <chrono>
#include <cstdio>
#include <iostream>

namespace evo_alg {
    template <class IndividualType, class FitnessType>
    std::tuple<IndividualType, Population<IndividualType>, std::vector<double>, std::vector<double>,
               std::vector<double>>
    ga(size_t const iterations, size_t const pop_size, size_t elite_size,
       typename FitnessType::const_shared_ptr const fitness,
       typename initializator::initialization_function_t<IndividualType, FitnessType> const init_fun,
       typename selector::selection_function_t<IndividualType> const selection_fun,
       typename recombinator::crossover_function_t<IndividualType> const crossover_fun, double const crossover_pr,
       typename mutator::mutation_function_t<IndividualType> const mutation_fun, double const mutation_pr,
       size_t const log_step = 0, double const convergence_threshold = NAN, size_t const convergence_interval = 100,
       double const convergence_eps = utils::eps) {
        std::vector<double> best_fit, mean_fit, diversity;

        std::chrono::duration<double, std::milli> it_time, gen_time, eval_time;
        std::chrono::time_point<std::chrono::high_resolution_clock> t1, t2, tt1, tt2;

        Population<IndividualType> population = init_fun(fitness, pop_size);
        population.evaluateFitness();

        best_fit.push_back(population.getBestFitness());
        mean_fit.push_back(population.getMeanFitness());
        diversity.push_back(population.getPairwiseDiversity());

        IndividualType best_individual = population[population.getBestIndividuals()[0]];
        size_t convergence_level = 0;
        double c = 1.2, c_inc = 0.8 / (double) iterations;
        for (size_t it = 0; it < iterations; ++it) {
            tt1 = std::chrono::high_resolution_clock::now();

            if (elite_size == 0 && it >= iterations / 2) {
                const std::vector<size_t> sorted_individuals = population.getSortedIndividuals();
                population[sorted_individuals[pop_size - 1]] = best_individual;
                elite_size = 1;
            }

            std::vector<double> pop_fitness(pop_size);
            for (size_t index = 0; index < pop_size; ++index)
                pop_fitness[index] = population.getIndividual(index).getFitnessValue()[0];

            pop_fitness =
                fitness::linearNormalization(pop_fitness, *std::min_element(pop_fitness.begin(), pop_fitness.end()),
                                             *std::max_element(pop_fitness.begin(), pop_fitness.end()));
            pop_fitness = fitness::linearScale(pop_fitness, c);

            Population<IndividualType> new_population(pop_size);

            t1 = std::chrono::high_resolution_clock::now();

            // clang-format off

            #pragma omp parallel for schedule(dynamic)

            //clang-format on
            for (size_t index = 0; index < pop_size; index += 2) {
                size_t const parent_1 = selection_fun(pop_fitness);

                std::vector<double> pop_fitness_2;
                std::vector<size_t> pop_indexes;
                for (size_t index_2 = 0; index_2 < pop_size; ++index_2) {
                    if (index_2 != parent_1) {
                        pop_fitness_2.push_back(pop_fitness[index_2]);
                        pop_indexes.push_back(index_2);
                    }
                }
                size_t const parent_2 = pop_indexes[selection_fun(pop_fitness_2)];

                double const pr = utils::uniformProbGen();
                IndividualType child_1;
                IndividualType child_2;
                if (pr < crossover_pr) {
                    std::tie(child_1, child_2) =
                        crossover_fun(population.getIndividual(parent_1), population.getIndividual(parent_2));
                } else {
                    child_1 = population.getIndividual(parent_1);
                    child_2 = population.getIndividual(parent_2);
                }

                child_1 = mutation_fun(child_1, mutation_pr);
                child_2 = mutation_fun(child_2, mutation_pr);

                new_population[index] = child_1;
                if (index + 1 < pop_size)
                    new_population[index + 1] = child_2;
            }
            t2 = std::chrono::high_resolution_clock::now();
            gen_time = t2 - t1;

            t1 = std::chrono::high_resolution_clock::now();
            new_population.evaluateFitness();
            t2 = std::chrono::high_resolution_clock::now();
            eval_time = t2 - t1;

            if (elite_size > 0) {
                const std::vector<size_t> sorted_individuals = population.getSortedIndividuals();
                const std::vector<size_t> new_sorted_individuals = new_population.getSortedIndividuals();
                for (size_t index = 0; index < elite_size; ++index)
                    new_population[new_sorted_individuals[pop_size - index - 1]] =
                        population[sorted_individuals[index]];
            }

            population = std::move(new_population);

            IndividualType current_best_individual = population[population.getBestIndividuals()[0]];
            if (utils::numericEqual(current_best_individual.getFitnessValue()[0], best_individual.getFitnessValue()[0],
                                    convergence_eps)) {
                ++convergence_level;
            } else {
                convergence_level = 0;
            }

            if (current_best_individual > best_individual)
                best_individual = current_best_individual;

            tt2 = std::chrono::high_resolution_clock::now();
            it_time = tt2 - tt1;

            c = c + c_inc;

            if (log_step > 0 && it % log_step == 0) {
                printf("Iteration %lu -> cur. fitness: %.9f | best fitness: %.9f | c: %.5f", it,
                       current_best_individual.getFitnessValue()[0], best_individual.getFitnessValue()[0], c);
                printf(" || times -> total: %.0fms | gen: %.0fms | eval: %.0fms\n", it_time.count(), gen_time.count(),
                       eval_time.count());
            }

            best_fit.push_back(population.getBestFitness());
            mean_fit.push_back(population.getMeanFitness());
            diversity.push_back(population.getPairwiseDiversity());

            if (!std::isnan(convergence_threshold) && utils::numericGreater(best_individual.getFitnessValue()[0], convergence_threshold))
                break;

            if (convergence_level >= convergence_interval)
                break;
        }

        return {best_individual, population, best_fit, mean_fit, diversity};
    }
}

#endif