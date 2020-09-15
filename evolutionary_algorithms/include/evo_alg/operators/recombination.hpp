#ifndef GUARD_H_EVO_ALG_RECOMBINATION
#define GUARD_H_EVO_ALG_RECOMBINATION

#include "../../../include/evo_alg/core.hpp"

#include <utility>

namespace evo_alg {
    namespace recombinator {
        template <class IndividualType>
        using crossover_function_t =
            std::function<std::pair<IndividualType, IndividualType>(IndividualType const&, IndividualType const&)>;

        std::pair<real_individual_t, real_individual_t>
        sbx(real_individual_t const& parent_1, real_individual_t const& parent_2, uint32_t const n, double const pr);

        std::pair<real_individual_t, real_individual_t> arithmetic(real_individual_t const& parent_1,
                                                                   real_individual_t const& parent_2);

        template <class GeneType>
        std::pair<Individual<GeneType>, Individual<GeneType>> onePoint(Individual<GeneType> const& parent_1,
                                                                       Individual<GeneType> const& parent_2) {
            Individual<GeneType> child_1(parent_1), child_2(parent_2);

            std::vector<GeneType> parent_1_chromosome = parent_1.getChromosome();
            std::vector<GeneType> parent_2_chromosome = parent_2.getChromosome();
            std::vector<GeneType> child_1_chromosome = child_1.getChromosome();
            std::vector<GeneType> child_2_chromosome = child_2.getChromosome();

            size_t const chromosome_size = parent_1.getChromosome().size();
            std::uniform_int_distribution<size_t> point_dist(1, chromosome_size - 1);
            size_t const cut_point = point_dist(utils::rng);
            for (size_t index = cut_point; index < chromosome_size; ++index) {
                child_1_chromosome[index] = parent_2_chromosome[index];
                child_2_chromosome[index] = parent_1_chromosome[index];
            }
            child_1.setChromosome(child_1_chromosome);
            child_2.setChromosome(child_2_chromosome);

            return {child_1, child_2};
        }

        template <class GeneType>
        std::pair<Individual<GeneType>, Individual<GeneType>> twoPoints(Individual<GeneType> const& parent_1,
                                                                        Individual<GeneType> const& parent_2) {
            Individual<GeneType> child_1(parent_1), child_2(parent_2);

            std::vector<GeneType> parent_1_chromosome = parent_1.getChromosome();
            std::vector<GeneType> parent_2_chromosome = parent_2.getChromosome();
            std::vector<GeneType> child_1_chromosome = child_1.getChromosome();
            std::vector<GeneType> child_2_chromosome = child_2.getChromosome();

            size_t const chromosome_size = parent_1.getChromosome().size();
            std::vector<size_t> indexes(chromosome_size);
            std::iota(indexes.begin(), indexes.end(), 0);
            std::shuffle(indexes.begin(), indexes.end(), utils::rng);

            size_t start_index = indexes[0];
            size_t end_index = indexes[1];

            if (start_index > end_index) {
                std::swap(start_index, end_index);
            }
            if (start_index == 0 && end_index == chromosome_size - 1) {
                end_index = indexes[2];
            }

            for (size_t index = start_index; index <= end_index; ++index) {
                child_1_chromosome[index] = parent_2_chromosome[index];
                child_2_chromosome[index] = parent_1_chromosome[index];
            }
            child_1.setChromosome(child_1_chromosome);
            child_2.setChromosome(child_2_chromosome);

            return {child_1, child_2};
        }
    }
}

#endif