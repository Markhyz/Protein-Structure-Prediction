#ifndef GUARD_H_PHENOTYPE
#define GUARD_H_PHENOTYPE

#include "../../../commons/include/macros.hpp"
#include "fitness.hpp"
#include <optional>
#include <stdexcept>

namespace EvoAlg {
    class AbstractPhenotype {
      public:
        POINTER_ALIAS(AbstractPhenotype)

        virtual ~AbstractPhenotype() = 0;
    };

    AbstractPhenotype::~AbstractPhenotype(){};

    template <size_t FitnessSize>
    class Phenotype : public AbstractPhenotype {
      public:
        POINTER_ALIAS(Phenotype<FitnessSize>)

        class UndefinedFitnessException : public std::exception {
          public:
            virtual const char* what() const throw() {
                return "undefined fitness value";
            }
        };

        Phenotype(typename AbstractFitnessFunction<FitnessSize>::shared_ptr const& fitness);

        void evaluateFitness(AbstractGenotype const& genotype);
        typename AbstractFitnessFunction<FitnessSize>::fitness_t getFitnessValue() const;

      private:
        typename AbstractFitnessFunction<FitnessSize>::shared_ptr fitness_;
        std::optional<typename AbstractFitnessFunction<FitnessSize>::fitness_t> fitness_value_;
    };

    template <size_t FitnessSize>
    Phenotype<FitnessSize>::Phenotype(typename AbstractFitnessFunction<FitnessSize>::shared_ptr const& fitness)
        : fitness_(fitness){};

    template <size_t FitnessSize>
    void Phenotype<FitnessSize>::evaluateFitness(AbstractGenotype const& genotype) {
        fitness_value_ = (*fitness_)(genotype);
    }

    template <size_t FitnessSize>
    typename AbstractFitnessFunction<FitnessSize>::fitness_t Phenotype<FitnessSize>::getFitnessValue() const {
        if (!fitness_value_.has_value()) {
            throw UndefinedFitnessException();
        }

        return fitness_value_.value();
    }
}

#endif