#ifndef GUARD_H_EVO_ALG_POPULATION
#define GUARD_H_EVO_ALG_POPULATION

#include "../commons/macros.hpp"
#include "individual.hpp"

#include <algorithm>

namespace evo_alg {
    template <typename... ChromosomeTypes>
    class Population {
      public:
        Population() = delete;
    };

    template <typename... ChromosomeTypes>
    class Population<Individual<ChromosomeTypes...>> {
      public:
        POINTER_ALIAS(Population)

        Population();
        Population(std::vector<Individual<ChromosomeTypes...>> const& individuals);
        Population(size_t const pop_size);

        size_t getSize() const;
        Individual<ChromosomeTypes...> const& getIndividual(size_t const index) const;

        void setIndividual(size_t const index, Individual<ChromosomeTypes...> const& individual);

        void appendIndividual(Individual<ChromosomeTypes...> const& Individual);

        void evaluateFitness();

        void resize(size_t const new_size);

        Individual<ChromosomeTypes...>& operator[](size_t const index);
        Individual<ChromosomeTypes...> const& operator[](size_t const index) const;

      private:
        std::vector<Individual<ChromosomeTypes...>> population_;
    };

    template <typename... ChromosomeTypes>
    Population<Individual<ChromosomeTypes...>>::Population(){};

    template <typename... ChromosomeTypes>
    Population<Individual<ChromosomeTypes...>>::Population(
        std::vector<Individual<ChromosomeTypes...>> const& individuals)
        : population_{individuals} {};

    template <typename... ChromosomeTypes>
    Population<Individual<ChromosomeTypes...>>::Population(size_t const pop_size) : population_{pop_size} {};

    template <typename... ChromosomeTypes>
    size_t Population<Individual<ChromosomeTypes...>>::getSize() const {
        return population_.size();
    }

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...> const&
    Population<Individual<ChromosomeTypes...>>::getIndividual(size_t const index) const {
        return (*this)[index];
    }

    template <typename... ChromosomeTypes>
    void Population<Individual<ChromosomeTypes...>>::setIndividual(size_t const index,
                                                                   Individual<ChromosomeTypes...> const& individual) {
        (*this)[index] = individual;
    }

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...>& Population<Individual<ChromosomeTypes...>>::operator[](size_t const index) {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <typename... ChromosomeTypes>
    Individual<ChromosomeTypes...> const&
    Population<Individual<ChromosomeTypes...>>::operator[](size_t const index) const {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <typename... ChromosomeTypes>
    void Population<Individual<ChromosomeTypes...>>::evaluateFitness() {
        std::for_each(population_.begin(), population_.end(),
                      [](Individual<ChromosomeTypes...>& ind) { ind.evaluateFitness(); });
    }

    template <typename... ChromosomeTypes>
    void
    Population<Individual<ChromosomeTypes...>>::appendIndividual(Individual<ChromosomeTypes...> const& individual) {
        population_.push_back(individual);
    }

    template <typename... ChromosomeTypes>
    void Population<Individual<ChromosomeTypes...>>::resize(size_t const new_size) {
        population_.resize(new_size);
    }
}

#endif