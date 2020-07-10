#ifndef GUARD_H_EVO_ALG_POPULATION
#define GUARD_H_EVO_ALG_POPULATION

#include "../commons.hpp"
#include "individual.hpp"

#include <algorithm>

namespace EvoAlg {
    template <class IndividualType>
    class Population {
      public:
        POINTER_ALIAS(Population)

        Population();
        Population(std::vector<IndividualType> const& individuals);
        Population(size_t pop_size);

        size_t getSize() const;
        IndividualType const& getIndividual(size_t index) const;

        void setIndividual(size_t index, IndividualType const& individual);

        void appendIndividual(IndividualType const& Individual);

        void evaluateFitness();

        IndividualType& operator[](size_t index);
        IndividualType const& operator[](size_t index) const;

      private:
        std::vector<IndividualType> population_;
    };

    template <class IndividualType>
    Population<IndividualType>::Population(){};

    template <class IndividualType>
    Population<IndividualType>::Population(std::vector<IndividualType> const& individuals) : population_(individuals){};

    template <class IndividualType>
    Population<IndividualType>::Population(size_t pop_size) : population_(pop_size){};

    template <class IndividualType>
    size_t Population<IndividualType>::getSize() const {
        return population_.size();
    }

    template <class IndividualType>
    IndividualType const& Population<IndividualType>::getIndividual(size_t index) const {
        return (*this)[index];
    }

    template <class IndividualType>
    void Population<IndividualType>::setIndividual(size_t index, IndividualType const& individual) {
        (*this)[index] = individual;
    }

    template <class IndividualType>
    IndividualType& Population<IndividualType>::operator[](size_t index) {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <class IndividualType>
    IndividualType const& Population<IndividualType>::operator[](size_t index) const {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <class IndividualType>
    void Population<IndividualType>::evaluateFitness() {
        std::for_each(population_.begin(), population_.end(), [](IndividualType& ind) { ind.evaluateFitness(); });
    }

    template <class IndividualType>
    void Population<IndividualType>::appendIndividual(IndividualType const& individual) {
        population_.push_back(individual);
    }
}

#endif