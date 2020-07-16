#ifndef GUARD_H_EVO_ALG_PHENOTYPE
#define GUARD_H_EVO_ALG_PHENOTYPE

#include "../commons/macros.hpp"
#include "fitness.hpp"

#include <optional>
#include <stdexcept>

namespace evo_alg {
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
        Phenotype(typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr const& fitness);

        void evaluateFitness(Genotype<ChromosomeTypes...> const& genotype);

        typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr getFitnessFunction() const;
        typename FitnessFunction<ChromosomeTypes...>::fitness_t const& getFitnessValue() const;

        template <size_t ChromosomeIndex = 0>
        std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                              types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
        getBounds() const;

        template <size_t ChromosomeIndex = 0, typename... Args>
        void setBounds(Args&&... args);

      private:
        typename FitnessFunction<ChromosomeTypes...>::shared_ptr fitness_;
        std::optional<typename FitnessFunction<ChromosomeTypes...>::fitness_t> fitness_value_;
    };

    template <typename... ChromosomeTypes>
    Phenotype<ChromosomeTypes...>::Phenotype(){};

    template <typename... ChromosomeTypes>
    Phenotype<ChromosomeTypes...>::Phenotype(
        typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr const& fitness)
        : fitness_(fitness->clone()){};

    template <typename... ChromosomeTypes>
    void Phenotype<ChromosomeTypes...>::evaluateFitness(Genotype<ChromosomeTypes...> const& genotype) {
        fitness_value_ = (*fitness_)(genotype);
    }

    template <typename... ChromosomeTypes>
    typename FitnessFunction<ChromosomeTypes...>::const_shared_ptr
    Phenotype<ChromosomeTypes...>::getFitnessFunction() const {
        return fitness_;
    }

    template <typename... ChromosomeTypes>
    typename FitnessFunction<ChromosomeTypes...>::fitness_t const&
    Phenotype<ChromosomeTypes...>::getFitnessValue() const {
        if (!fitness_value_.has_value()) {
            throw UndefinedFitnessException();
        }

        return fitness_value_.value();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                          types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
    Phenotype<ChromosomeTypes...>::getBounds() const {
        return fitness_->template getBounds<ChromosomeIndex>();
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex, typename... Args>
    void Phenotype<ChromosomeTypes...>::setBounds(Args&&... args) {
        fitness_->template setBounds<ChromosomeIndex>(std::forward<Args>(args)...);
    }
}

#endif