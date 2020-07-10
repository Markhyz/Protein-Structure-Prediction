#ifndef GUARD_H_EVO_ALG_PHENOTYPE
#define GUARD_H_EVO_ALG_PHENOTYPE

#include "../commons.hpp"
#include "fitness.hpp"

#include <optional>
#include <stdexcept>

namespace EvoAlg {
    template <typename... ChromosomeTypes>
    class Phenotype {
      public:
        POINTER_ALIAS(Phenotype)

        class UndefinedFitnessException : public std::exception {
          public:
            virtual const char* what() const throw() {
                return "undefined fitness value";
            }
        };

        Phenotype();
        Phenotype(typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr const& fitness);

        void evaluateFitness(Genotype<ChromosomeTypes...> const& genotype);

        typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr getFitnessFunction() const;
        typename AbstractFitnessFunction<ChromosomeTypes...>::fitness_t const& getFitnessValue() const;

      private:
        typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr fitness_;
        std::optional<typename AbstractFitnessFunction<ChromosomeTypes...>::fitness_t> fitness_value_;
    };

    template <typename... ChromosomeTypes>
    Phenotype<ChromosomeTypes...>::Phenotype(){};

    template <typename... ChromosomeTypes>
    Phenotype<ChromosomeTypes...>::Phenotype(
        typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr const& fitness)
        : fitness_(fitness){};

    template <typename... ChromosomeTypes>
    void Phenotype<ChromosomeTypes...>::evaluateFitness(Genotype<ChromosomeTypes...> const& genotype) {
        fitness_value_ = (*fitness_)(genotype);
    }

    template <typename... ChromosomeTypes>
    typename AbstractFitnessFunction<ChromosomeTypes...>::const_shared_ptr
    Phenotype<ChromosomeTypes...>::getFitnessFunction() const {
        return fitness_;
    }

    template <typename... ChromosomeTypes>
    typename AbstractFitnessFunction<ChromosomeTypes...>::fitness_t const&
    Phenotype<ChromosomeTypes...>::getFitnessValue() const {
        if (!fitness_value_.has_value()) {
            throw UndefinedFitnessException();
        }

        return fitness_value_.value();
    }
}

#endif