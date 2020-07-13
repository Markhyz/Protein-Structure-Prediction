#ifndef GUARD_H_EVO_ALG_INDIVIDUAL
#define GUARD_H_EVO_ALG_INDIVIDUAL

#include "../commons/macros.hpp"
#include "../commons/types.hpp"
#include "genotype.hpp"
#include "phenotype.hpp"

#include <utility>

namespace evo_alg {
    template <typename... ChromosomeTypes>
    class Individual {
      public:
        POINTER_ALIAS(Individual)

        Individual();
        Individual(typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr fitness,
                   std::vector<ChromosomeTypes> const&... chromosomes);

        void evaluateFitness();

        template <size_t ChromosomeIndex = 0>
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> getChromosome() const;

        typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr getFitnessFunction() const;
        typename AbstractFitnessFunction<ChromosomeTypes...>::fitness_t const& getFitnessValue() const;

        template <size_t ChromosomeIndex, typename... Args>
        void setChromosome(Args&&... args);

      private:
        Genotype<ChromosomeTypes...> genotype_;
        Phenotype<ChromosomeTypes...> phenotype_;
    };

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...>::Individual(){};

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...>::Individual(
        typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr fitness,
        std::vector<ChromosomeTypes> const&... chromosomes)
        : genotype_(chromosomes...), phenotype_(fitness){};

    template <typename... ChromosomeTypes>
    void Individual<ChromosomeTypes...>::evaluateFitness() {
        phenotype_.evaluateFitness(genotype_);
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> Individual<ChromosomeTypes...>::getChromosome() const {
        return genotype_.template getChromosome<ChromosomeIndex>();
    }

    template <typename... ChromosomeTypes>
    typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr
    Individual<ChromosomeTypes...>::getFitnessFunction() const {
        return phenotype_.getFitnessFunction();
    }

    template <typename... ChromosomeTypes>
    typename AbstractFitnessFunction<ChromosomeTypes...>::fitness_t const&
    Individual<ChromosomeTypes...>::getFitnessValue() const {
        return phenotype_.getFitnessValue();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex, typename... Args>
    void Individual<ChromosomeTypes...>::setChromosome(Args&&... args) {
        genotype_.template setChromosome<ChromosomeIndex>(std::forward<Args>(args)...);
    }
}

#endif