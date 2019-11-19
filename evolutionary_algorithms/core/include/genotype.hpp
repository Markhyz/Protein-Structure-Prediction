#ifndef GUARD_H_GENOTYPE
#define GUARD_H_GENOTYPE

#include "../../../commons/include/macros.hpp"
#include "../../../commons/include/types.hpp"
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

constexpr bool cmp(int x, int y) {
    return x < y;
}

namespace EvoAlg {
    class AbstractGenotype {
      public:
        virtual ~AbstractGenotype() = 0;
    };

    AbstractGenotype::~AbstractGenotype(){};

    template <typename... ChromosomeTypes>
    class Genotype : public AbstractGenotype {
      public:
        POINTER_ALIAS(Genotype)

        template <size_t ChromosomeIndex>
        using ValidChromosomeIndex = std::enable_if_t<(ChromosomeIndex < sizeof...(ChromosomeTypes)), int>;

        using gene_t = std::variant<ChromosomeTypes...>;
        using chromosome_t = std::vector<gene_t>;

        template <size_t ChromosomeIndex>
        using NthGeneType = NthType<ChromosomeIndex, ChromosomeTypes...>;

        Genotype(std::vector<ChromosomeTypes> const&... chromosomes);
        Genotype(std::vector<size_t> const& chromosome_size);

        template <size_t ChromosomeIndex, ValidChromosomeIndex<ChromosomeIndex> = 0>
        void setChromosomeGene(size_t index, NthGeneType<ChromosomeIndex> const& gene);

        template <size_t ChromosomeIndex, ValidChromosomeIndex<ChromosomeIndex> = 0>
        NthGeneType<ChromosomeIndex> getChromosomeGene(size_t index) const;

        template <size_t ChromosomeIndex, ValidChromosomeIndex<ChromosomeIndex> = 0>
        std::vector<NthGeneType<ChromosomeIndex>> getChromosome() const;

      private:
        template <typename ChromosomeType>
        chromosome_t createChromosome(std::vector<ChromosomeType> const& chromosome) const;

        std::vector<chromosome_t> chromosomes_;
        std::vector<size_t> chromosome_size_;
    };

    template <typename... ChromosomeTypes>
    Genotype<ChromosomeTypes...>::Genotype(std::vector<ChromosomeTypes> const&... chromosomes)
        : chromosomes_{this->createChromosome(chromosomes)...}, chromosome_size_(chromosomes_.size()) {
        for (size_t index = 0; index < chromosomes_.size(); ++index) {
            chromosome_size_[index] = chromosomes_[index].size();
        }
    }

    template <typename... ChromosomeTypes>
    Genotype<ChromosomeTypes...>::Genotype(std::vector<size_t> const& chromosome_size)
        : chromosomes_(sizeof...(ChromosomeTypes)), chromosome_size_(chromosome_size) {
        if (chromosome_size_.size() > sizeof...(ChromosomeTypes)) {
            throw std::length_error("number of sizes greater than number of types");
        }
        for (size_t index = 0; index < chromosomes_.size(); ++index) {
            chromosomes_[index].resize(chromosome_size_[index]);
        }
    }

    template <typename... ChromosomeTypes>
    template <typename ChromosomeType>
    typename Genotype<ChromosomeTypes...>::chromosome_t
    Genotype<ChromosomeTypes...>::createChromosome(std::vector<ChromosomeType> const& chromosome) const {
        chromosome_t variant_chromosome(chromosome.size());
        std::copy(chromosome.begin(), chromosome.end(), variant_chromosome.begin());
        return variant_chromosome;
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex,
              typename Genotype<ChromosomeTypes...>::template ValidChromosomeIndex<ChromosomeIndex>>
    void Genotype<ChromosomeTypes...>::setChromosomeGene(
        size_t index, typename Genotype<ChromosomeTypes...>::template NthGeneType<ChromosomeIndex> const& gene) {
        chromosomes_[ChromosomeIndex][index] = gene;
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex,
              typename Genotype<ChromosomeTypes...>::template ValidChromosomeIndex<ChromosomeIndex>>
    typename Genotype<ChromosomeTypes...>::template NthGeneType<ChromosomeIndex>
    Genotype<ChromosomeTypes...>::getChromosomeGene(size_t index) const {
        return chromosomes_[ChromosomeIndex][index];
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex,
              typename Genotype<ChromosomeTypes...>::template ValidChromosomeIndex<ChromosomeIndex>>
    std::vector<typename Genotype<ChromosomeTypes...>::template NthGeneType<ChromosomeIndex>>
    Genotype<ChromosomeTypes...>::getChromosome() const {
        chromosome_t chromosome = chromosomes_[ChromosomeIndex];
        std::vector<NthGeneType<ChromosomeIndex>> simple_chromosome(chromosome.size());
        transform(chromosome.begin(), chromosome.end(), simple_chromosome.begin(),
                  [](gene_t const& gene) -> NthGeneType<ChromosomeIndex> { return std::get<ChromosomeIndex>(gene); });
        return simple_chromosome;
    }
}

#endif