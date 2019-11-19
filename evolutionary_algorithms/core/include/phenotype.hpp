#ifndef GUARD_H_PHENOTYPE
#define GUARD_H_PHENOTYPE

#include "../../../commons/include/macros.hpp"
#include "fitness.hpp"
#include <optional>
#include <stdexcept>

namespace EvoAlg {
    class Phenotype {
      public:
        POINTER_ALIAS(Phenotype)

        class UndefinedFitnessException : public std::exception {
          public:
            virtual const char* what() const throw();
        };

        Phenotype(AbstractFitness::shared_ptr const& fitness);

        double getFitnessValue(size_t index) const;
        void evaluateFitness(AbstractGenotype const& genotype);

      private:
        AbstractFitness::shared_ptr fitness_;
        std::optional<AbstractFitness::fitness_t> fitness_value_;
    };
}

#endif