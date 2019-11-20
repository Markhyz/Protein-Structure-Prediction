#ifndef GUARD_H_INDIVIDUAL
#define GUARD_H_INDIVIDUAL

#include "../../../commons/include/macros.hpp"
#include "../../../commons/include/types.hpp"
#include "genotype.hpp"
#include "phenotype.hpp"

namespace EvoAlg {
    class AbstractIndividual {
      public:
        POINTER_ALIAS(AbstractIndividual)

        virtual ~AbstractIndividual() = 0;
    };

    AbstractIndividual::~AbstractIndividual(){};

    template <typename... ChromosomeTypes>
    class Individual : public AbstractIndividual {
      public:
        POINTER_ALIAS(Individual)

        Individual(AbstractFitness::shared_ptr fitness, std::vector<ChromosomeTypes> const&... chromosomes);
        Individual(AbstractFitness::shared_ptr fitness, std::vector<size_t> const& chromosome_size);

        void evaluateFitness();

        template <size_t ChromosomeIndex>
        std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> getChromosome() const;

        double getFitnessValue(size_t index) const;

        template <size_t ChromosomeIndex>
        void setChromosome(std::vector<NthType<ChromosomeIndex, ChromosomeTypes...>> const& chromosome);

        template <size_t ChromosomeIndex>
        void setChromosome(size_t index, NthType<ChromosomeIndex, ChromosomeTypes...> const& gene);

      private:
        Genotype<ChromosomeTypes...> genotype;
        Phenotype phenotype;
    };
}

#endif