#ifndef GUARD_H_EVO_ALG_GENOTYPE
#define GUARD_H_EVO_ALG_GENOTYPE

#include "../commons/macros.hpp"
#include "../commons/types.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace evo_alg {
    template <typename... ChromosomeTypes>
    class Genotype {
      public:
        POINTER_ALIAS(Genotype<ChromosomeTypes...>)

        using gene_t = std::variant<ChromosomeTypes...>;
        using chromosome_t = std::vector<gene_t>;

        Genotype();
        Genotype(std::vector<ChromosomeTypes> const&... chromosomes);

        template <size_t ChromosomeIndex = 0>
        std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>> getChromosome() const;

        template <size_t ChromosomeIndex = 0>
        void setChromosome(std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>> const& chromosome);

        template <size_t ChromosomeIndex = 0>
        void setChromosome(size_t const index, types::NthType<ChromosomeIndex, ChromosomeTypes...> const& gene);

      private:
        std::vector<chromosome_t> chromosomes_;
    };

    template <typename... ChromosomeTypes>
    Genotype<ChromosomeTypes...>::Genotype(){};

    template <typename... ChromosomeTypes>
    Genotype<ChromosomeTypes...>::Genotype(std::vector<ChromosomeTypes> const&... chromosomes)
        : chromosomes_{
              typename Genotype<ChromosomeTypes...>::chromosome_t(chromosomes.begin(), chromosomes.end())...} {};

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>>
    Genotype<ChromosomeTypes...>::getChromosome() const {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>> chromosome;
        std::transform(chromosomes_[ChromosomeIndex].begin(), chromosomes_[ChromosomeIndex].end(),
                       std::back_inserter(chromosome),
                       [](gene_t const& gene) { return std::get<ChromosomeIndex>(gene); });

        return chromosome;
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void Genotype<ChromosomeTypes...>::setChromosome(
        std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>> const& chromosome) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        chromosomes_[ChromosomeIndex] =
            typename Genotype<ChromosomeTypes...>::chromosome_t(chromosome.begin(), chromosome.end());
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void Genotype<ChromosomeTypes...>::setChromosome(size_t const index,
                                                     types::NthType<ChromosomeIndex, ChromosomeTypes...> const& gene) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        if (index >= chromosomes_[ChromosomeIndex].size()) {
            throw std::out_of_range("out of range access");
        }

        chromosomes_[ChromosomeIndex][index] = gene;
    }
}

#endif