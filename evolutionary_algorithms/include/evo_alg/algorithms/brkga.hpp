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
    namespace brkga {
        using decoding_function_t = std::function<std::vector<double>(std::vector<double> const&)>;

        template <class FitnessType>
        class BrkgaFitness : public FitnessFunction<double> {
          public:
            BrkgaFitness(size_t chromosome_size, typename FitnessType::const_shared_ptr fitness,
                         decoding_function_t decoder)
                : FitnessFunction<double>({chromosome_size, {0, 1}}), chromosome_size_(chromosome_size),
                  fitness_(fitness->clone()), decoder_(decoder){};

            size_t getDimension() const override {
                return dimension_;
            }

            BrkgaFitness* clone() const override {
                return new BrkgaFitness(chromosome_size_, fitness_, decoder_);
            }

            fitness_t operator()(evo_alg::Genotype<double> const& genotype) override {
                std::vector<double> encoded_chromosome = genotype.getChromosome();

                std::vector<double> result = (*fitness_)({decoder_(encoded_chromosome)});

                return result;
            }

          private:
            size_t chromosome_size_;
            typename FitnessType::shared_ptr fitness_;
            decoding_function_t decoder_;
            size_t dimension_ = 1;
        };

        template <class FitnessType>
        struct config_t {
            config_t(size_t iteration_num, size_t pop_size, size_t chromosome_size, double elite_fraction,
                     double mut_fraction, double elite_cross_pr, typename FitnessType::const_shared_ptr fitness,
                     decoding_function_t decoder)
                : iteration_num(iteration_num), pop_size(pop_size), chromosome_size(chromosome_size),
                  elite_fraction(elite_fraction), mut_fraction(mut_fraction), elite_cross_pr(elite_cross_pr),
                  fitness(fitness), decoder(decoder){};

            size_t iteration_num;
            size_t pop_size;
            size_t chromosome_size;
            double elite_fraction;
            double mut_fraction;
            double elite_cross_pr;
            typename FitnessType::const_shared_ptr fitness;
            decoding_function_t decoder;
            size_t log_step = 0;
            double diversity_threshold = 0;
            double convergence_threshold = NAN;
            size_t convergence_interval = 0;
            double convergence_eps = utils::eps;
            std::vector<std::vector<double>> initial_pop;
            std::function<void(config_t&, size_t)> update_fn;
        };

        template <class IndividualType, class FitnessType>
        std::tuple<IndividualType, Population<IndividualType>, std::vector<double>, std::vector<double>,
                   std::vector<double>>
        run(config_t<FitnessType> config) {
            size_t const iteration_num = config.iteration_num;
            size_t const pop_size = config.pop_size;
            size_t const chromosome_size = config.chromosome_size;
            typename FitnessType::const_shared_ptr const fitness = config.fitness;
            decoding_function_t decoder = config.decoder;
            std::vector<std::vector<double>> const& initial_pop = config.initial_pop;

            FitnessFunction<double>::shared_ptr brkga_fitness(
                new BrkgaFitness<FitnessType>(chromosome_size, fitness, decoder));

            size_t elite_size = (size_t)((double) pop_size * config.elite_fraction);
            size_t mut_size = (size_t)((double) pop_size * config.mut_fraction);
            size_t cross_size = pop_size - elite_size - mut_size;

            std::vector<double> best_fit, mean_fit, diversity;

            std::chrono::duration<double, std::milli> it_time, gen_time, eval_time;
            std::chrono::time_point<std::chrono::high_resolution_clock> t1, t2, tt1, tt2;

            Population<real_individual_t> population(pop_size), new_population(pop_size);
            Population<IndividualType> decoded_population(pop_size);

#pragma omp parallel for schedule(dynamic)
            for (size_t index = 0; index < pop_size; ++index) {
                population[index] = {brkga_fitness};
                new_population[index] = {brkga_fitness};
                decoded_population[index] = {fitness};
            }
            if (pop_size - initial_pop.size() > 0) {
                initializator::uniformRandomInit(population, brkga_fitness, pop_size - initial_pop.size());
            }
            for (size_t index = 0; index < initial_pop.size(); ++index) {
                population[pop_size - initial_pop.size() + index].setChromosome(initial_pop[index]);
            }

            population.evaluateFitness();

            real_individual_t child(population[0]);

            best_fit.push_back(population.getBestFitness());
            mean_fit.push_back(population.getMeanFitness());
            diversity.push_back(population.getPairwiseDiversity());

            real_individual_t best_individual = population[population.getBestIndividuals()[0]];
            size_t convergence_level = 0;
            for (size_t it = 0; it < iteration_num; ++it) {
                tt1 = std::chrono::high_resolution_clock::now();

                std::vector<double> pop_fitness(pop_size);
                for (size_t index = 0; index < pop_size; ++index)
                    pop_fitness[index] = population[index].getFitnessValue()[0];

                initializator::uniformRandomInit(new_population, brkga_fitness, mut_size);

                std::vector<size_t> sorted_individuals = population.getSortedIndividuals();
                std::vector<size_t> diversified_individuals, remaining_individuals;
                for (size_t index = 0; index < sorted_individuals.size(); ++index) {
                    if (diversified_individuals.size() < elite_size) {
                        real_individual_t const& ind_1 = population[sorted_individuals[index]];
                        size_t index_2 = 0;
                        for (; index_2 < diversified_individuals.size(); ++index_2) {
                            real_individual_t const& ind_2 = population[diversified_individuals[index_2]];
                            double const distance = ind_1.getEuclidianDistance(ind_2);
                            if (distance < config.diversity_threshold)
                                break;
                        }
                        if (index_2 == diversified_individuals.size()) {
                            diversified_individuals.push_back(sorted_individuals[index]);
                            continue;
                        }
                    }
                    remaining_individuals.push_back(sorted_individuals[index]);
                }
                for (size_t index = 0; index < remaining_individuals.size(); ++index) {
                    diversified_individuals.push_back(remaining_individuals[index]);
                }

                for (size_t index = 0; index < elite_size; ++index) {
                    const size_t ind_index = diversified_individuals[index];
                    new_population[mut_size + index].setChromosome(population[ind_index].getChromosome());
                }

                t1 = std::chrono::high_resolution_clock::now();
                for (size_t index = 0; index < cross_size; ++index) {
                    std::uniform_int_distribution<size_t> elite_dist(0, elite_size - 1);
                    std::uniform_int_distribution<size_t> non_elite_dist(elite_size, pop_size - 1);

                    size_t const elite_parent = diversified_individuals[elite_dist(utils::rng)];
                    size_t const non_elite_parent = diversified_individuals[non_elite_dist(utils::rng)];

                    recombinator::eliteUniform<double>(population[elite_parent], population[non_elite_parent], child,
                                                       config.elite_cross_pr);

                    new_population[mut_size + elite_size + index].setChromosome(child.getChromosome());
                }
                t2 = std::chrono::high_resolution_clock::now();
                gen_time = t2 - t1;

                assert(mut_size + elite_size + cross_size == pop_size);

                t1 = std::chrono::high_resolution_clock::now();
                new_population.evaluateFitness();
                t2 = std::chrono::high_resolution_clock::now();
                eval_time = t2 - t1;

                for (size_t index = 0; index < pop_size; ++index) {
                    population[index].setChromosome(new_population[index].getChromosome());
                    population[index].setFitnessValue(new_population[index].getFitnessValue());
                }

                size_t const current_best_individual = population.getBestIndividuals()[0];
                if (utils::numericEqual(population[current_best_individual].getFitnessValue()[0],
                                        best_individual.getFitnessValue()[0], config.convergence_eps)) {
                    ++convergence_level;
                } else {
                    convergence_level = 0;
                }

                if (population[current_best_individual] > best_individual) {
                    best_individual.setChromosome(population[current_best_individual].getChromosome());
                    best_individual.setFitnessValue(population[current_best_individual].getFitnessValue());
                }

                double const best_fitness = best_individual.getFitnessValue()[0];
                best_fit.push_back(best_fitness);

                double const mean_fitness = population.getMeanFitness(mut_size);
                mean_fit.push_back(mean_fitness);

                double const diver = population.getPairwiseDiversity(mut_size);
                diversity.push_back(diver);

                tt2 = std::chrono::high_resolution_clock::now();
                it_time = tt2 - tt1;

                if (config.log_step > 0 && it % config.log_step == 0) {
                    printf(
                        "Gen %lu -> best fitness: %.7f | mean fitness: %.7f | diversity: %.3f | e: %.2f | m: %.2f | c: "
                        "%.2f",
                        it, best_fitness, mean_fitness, diver, config.elite_fraction, config.mut_fraction,
                        config.elite_cross_pr);
                    printf(" || it: %.0fms | gen: %.0fms | eval: %.0fms\n", it_time.count(), gen_time.count(),
                           eval_time.count());
                }

                if (config.update_fn)
                    config.update_fn(config, it);

                elite_size = (size_t)((double) pop_size * config.elite_fraction);
                mut_size = (size_t)((double) pop_size * config.mut_fraction);
                cross_size = pop_size - elite_size - mut_size;

                if (!std::isnan(config.convergence_threshold) &&
                    utils::numericGreater(best_individual.getFitnessValue()[0], config.convergence_threshold))
                    break;

                if (config.convergence_interval && convergence_level >= config.convergence_interval)
                    break;
            }

            IndividualType decoded_best_individual = {fitness, decoder(best_individual.getChromosome())};
            decoded_best_individual.setFitnessValue(best_individual.getFitnessValue());

            for (size_t index = 0; index < pop_size; ++index) {
                decoded_population[index].setChromosome(decoder(population[index].getChromosome()));
                decoded_population[index].setFitnessValue(population[index].getFitnessValue());
            }

            return {decoded_best_individual, decoded_population, best_fit, mean_fit, diversity};
        }
    }
}

#endif