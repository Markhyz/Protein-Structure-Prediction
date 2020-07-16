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
        Individual(typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr fitness,
                   std::vector<ChromosomeTypes> const&... chromosomes);

        void evaluateFitness();

        template <size_t ChromosomeIndex = 0>
        std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>> getChromosome() const;

        typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr getFitnessFunction() const;
        typename FitnessFunction<ChromosomeTypes...>::fitness_t const& getFitnessValue() const;

        template <size_t ChromosomeIndex = 0>
        std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                              types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
        getBounds() const;

        template <size_t ChromosomeIndex = 0, typename... Args>
        void setBounds(Args&&... args);

        template <size_t ChromosomeIndex = 0, typename... Args>
        void setChromosome(Args&&... args);

      private:
        Genotype<ChromosomeTypes...> genotype_;
        Phenotype<ChromosomeTypes...> phenotype_;
    };

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...>::Individual(){};

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...>::Individual(typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr fitness,
                                               std::vector<ChromosomeTypes> const&... chromosomes)
        : genotype_(chromosomes...), phenotype_(fitness){};

    template <typename... ChromosomeTypes>
    void Individual<ChromosomeTypes...>::evaluateFitness() {
        phenotype_.evaluateFitness(genotype_);
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<types::NthType<ChromosomeIndex, ChromosomeTypes...>>
    Individual<ChromosomeTypes...>::getChromosome() const {
        return genotype_.template getChromosome<ChromosomeIndex>();
    }

    template <typename... ChromosomeTypes>
    typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr
    Individual<ChromosomeTypes...>::getFitnessFunction() const {
        return phenotype_.getFitnessFunction();
    }

    template <typename... ChromosomeTypes>
    typename FitnessFunction<ChromosomeTypes...>::fitness_t const&
    Individual<ChromosomeTypes...>::getFitnessValue() const {
        return phenotype_.getFitnessValue();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                          types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
    Individual<ChromosomeTypes...>::getBounds() const {
        return phenotype_.template getBounds<ChromosomeIndex>();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex, typename... Args>
    void Individual<ChromosomeTypes...>::setBounds(Args&&... args) {
        phenotype_.template setBounds<ChromosomeIndex>(std::forward<Args>(args)...);
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex, typename... Args>
    void Individual<ChromosomeTypes...>::setChromosome(Args&&... args) {
        genotype_.template setChromosome<ChromosomeIndex>(std::forward<Args>(args)...);
    }
}

#endif