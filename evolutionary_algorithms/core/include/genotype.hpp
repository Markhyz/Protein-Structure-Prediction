#ifndef GUARD_H_GENOTYPE
#define GUARD_H_GENOTYPE

#include "../../commons/include/macros.hpp"
#include "../../commons/include/types.hpp"
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace EvoAlg {
    class AbstractGenotype {
      public:
        POINTER_ALIAS(AbstractGenotype)

        virtual ~AbstractGenotype() = 0;
    };

    template <typename... ChromosomeTypes>
    class Genotype : public AbstractGenotype {
      public:
        POINTER_ALIAS(Genotype<ChromosomeTypes...>)

        using gene_t = std::variant<ChromosomeTypes...>;
        using chromosome_t = std::vector<gene_t>;

        Genotype(std::vector<ChromosomeTypes> const&... chromosomes);

        template <size_t ChromosomeIndex>
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> getChromosome() const;

        template <size_t ChromosomeIndex>
        void setChromosome(std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> const& chromosome);

        template <size_t ChromosomeIndex>
        void setChromosome(size_t index, NthType<ChromosomeIndex, ChromosomeTypes...> const& gene);

      private:
        template <size_t ChromosomeIndex>
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>>
        decodeChromosome(std::vector<gene_t> const& encoded_chromosome) const;

        template <typename ChromosomeType>
        chromosome_t encodeChromosome(std::vector<ChromosomeType> const& simple_chromosome) const;

        std::vector<chromosome_t> chromosomes_;
        std::vector<size_t> chromosome_size_;
    };

    template <typename... ChromosomeTypes>
    Genotype<ChromosomeTypes...>::Genotype(std::vector<ChromosomeTypes> const&... chromosomes)
        : chromosomes_{this->encodeChromosome(chromosomes)...}, chromosome_size_(chromosomes_.size()) {
        for (size_t index = 0; index < chromosomes_.size(); ++index) {
            chromosome_size_[index] = chromosomes_[index].size();
        }
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> Genotype<ChromosomeTypes...>::decodeChromosome(
        std::vector<typename Genotype<ChromosomeTypes...>::gene_t> const& encoded_chromosome) const {
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> simple_chromosome(encoded_chromosome.size());
        transform(
            encoded_chromosome.begin(), encoded_chromosome.end(), simple_chromosome.begin(),
            [](gene_t const& gene) -> NthType<ChromosomeIndex, ChromosomeTypes...> {
                return std::get<ChromosomeIndex>(gene);
            });
        return simple_chromosome;
    }

    template <typename... ChromosomeTypes>
    template <typename ChromosomeType>
    typename Genotype<ChromosomeTypes...>::chromosome_t
    Genotype<ChromosomeTypes...>::encodeChromosome(std::vector<ChromosomeType> const& simple_chromosome) const {
        chromosome_t encoded_chromosome(simple_chromosome.size());
        std::copy(simple_chromosome.begin(), simple_chromosome.end(), encoded_chromosome.begin());
        return encoded_chromosome;
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> Genotype<ChromosomeTypes...>::getChromosome() const {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        chromosome_t chromosome = chromosomes_[ChromosomeIndex];
        return this->decodeChromosome<ChromosomeIndex>(chromosome);
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void Genotype<ChromosomeTypes...>::setChromosome(
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> const& chromosome) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        chromosomes_[ChromosomeIndex] = this->encodeChromosome(chromosome);
        chromosome_size_[ChromosomeIndex] = chromosomes_[ChromosomeIndex].size();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void Genotype<ChromosomeTypes...>::setChromosome(size_t index,
                                                     NthType<ChromosomeIndex, ChromosomeTypes...> const& gene) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        if (index >= chromosome_size_[ChromosomeIndex]) {
            throw std::out_of_range("out of range access");
        }

        chromosomes_[ChromosomeIndex][index] = gene;
    }
}

#endif