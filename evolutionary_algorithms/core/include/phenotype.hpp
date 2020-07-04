#ifndef GUARD_H_PHENOTYPE
#define GUARD_H_PHENOTYPE

#include "../../commons/include/macros.hpp"
#include "fitness.hpp"
#include <optional>
#include <stdexcept>

namespace EvoAlg {
    class Phenotype {
      public:
        POINTER_ALIAS(Phenotype)

        class UndefinedFitnessException : public std::exception {
          public:
            virtual const char* what() const throw() {
                return "undefined fitness value";
            }
        };

        Phenotype(typename AbstractFitnessFunction::const_shared_ptr const& fitness);

        void evaluateFitness(AbstractGenotype const& genotype);

        AbstractFitnessFunction::const_shared_ptr getFitnessFunction() const;
        typename AbstractFitnessFunction::fitness_t const& getFitnessValue() const;

      private:
        typename AbstractFitnessFunction::const_shared_ptr fitness_;
        std::optional<typename AbstractFitnessFunction::fitness_t> fitness_value_;
    };
}

#endif