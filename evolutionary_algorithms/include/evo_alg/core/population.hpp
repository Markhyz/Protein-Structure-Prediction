#ifndef GUARD_H_EVO_ALG_POPULATION
#define GUARD_H_EVO_ALG_POPULATION

#include "../commons/macros.hpp"
#include "../commons/types.hpp"
#include "individual.hpp"

#include <omp.h>

#include <algorithm>
#include <iostream>

namespace evo_alg {
    template <class IndividualType>
    class Population {
      public:
        POINTER_ALIAS(Population)

        Population();
        Population(std::vector<IndividualType> const individuals);
        Population(size_t const pop_size);

        size_t getSize() const;
        IndividualType const& getIndividual(size_t const index) const;
        std::vector<size_t> getBestIndividuals() const;
        std::vector<size_t> getSortedIndividuals() const;

        void setIndividual(size_t const index, IndividualType const& individual);

        void appendIndividual(IndividualType const& Individual);

        void evaluateFitness();

        void resize(size_t const new_size);

        double getBestFitness() const;
        double getMeanFitness() const;
        double getPairwiseDiversity() const;

        IndividualType& operator[](size_t const index);
        IndividualType const& operator[](size_t const index) const;

      protected:
        std::vector<IndividualType> population_;
    };

    template <class IndividualType>
    Population<IndividualType>::Population(){};

    template <class IndividualType>
    Population<IndividualType>::Population(std::vector<IndividualType> const individuals) : population_{individuals} {};

    template <class IndividualType>
    Population<IndividualType>::Population(size_t const pop_size) : population_{pop_size} {};

    template <class IndividualType>
    size_t Population<IndividualType>::getSize() const {
        return population_.size();
    }

    template <class IndividualType>
    IndividualType const& Population<IndividualType>::getIndividual(size_t const index) const {
        return (*this)[index];
    }

    template <class IndividualType>
    std::vector<size_t> Population<IndividualType>::getBestIndividuals() const {
        std::vector<size_t> best_individuals;
        size_t best_individual = 0;
        for (size_t index = 1; index < population_.size(); ++index)
            if (population_[index] > population_[best_individual])
                best_individual = index;
        for (size_t index = 0; index < population_.size(); ++index)
            if (!(population_[index] < population_[best_individual]))
                best_individuals.push_back(index);

        return best_individuals;
    }

    template <class IndividualType>
    std::vector<size_t> Population<IndividualType>::getSortedIndividuals() const {
        std::vector<size_t> sorted_individuals(population_.size());
        std::iota(sorted_individuals.begin(), sorted_individuals.end(), 0);
        std::sort(sorted_individuals.rbegin(), sorted_individuals.rend(),
                  [this](size_t const ind_1, size_t const ind_2) { return population_[ind_1] < population_[ind_2]; });

        return sorted_individuals;
    }

    template <class IndividualType>
    void Population<IndividualType>::setIndividual(size_t const index, IndividualType const& individual) {
        (*this)[index] = individual;
    }

    template <class IndividualType>
    IndividualType& Population<IndividualType>::operator[](size_t const index) {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <class IndividualType>
    IndividualType const& Population<IndividualType>::operator[](size_t const index) const {
        if (index >= population_.size()) {
            throw std::out_of_range("out of range access");
        }

        return population_[index];
    }

    template <class IndividualType>
    void Population<IndividualType>::evaluateFitness() {
        // #pragma omp parallel for schedule(dynamic)
        for (size_t index = 0; index < population_.size(); ++index)
            population_[index].evaluateFitness();
    }

    template <class IndividualType>
    void Population<IndividualType>::appendIndividual(IndividualType const& individual) {
        population_.push_back(individual);
    }

    template <class IndividualType>
    void Population<IndividualType>::resize(size_t const new_size) {
        population_.resize(new_size);
    }

    template <class IndividualType>
    double Population<IndividualType>::getBestFitness() const {
        double best_fit = this->population_[0].getFitnessValue()[0];
        for (size_t index = 1; index < this->population_.size(); ++index)
            best_fit = std::max(best_fit, this->population_[index].getFitnessValue()[0]);

        return best_fit;
    }

    template <class IndividualType>
    double Population<IndividualType>::getMeanFitness() const {
        double total_fit = 0;
        for (size_t index = 0; index < this->population_.size(); ++index)
            total_fit += this->population_[index].getFitnessValue()[0];

        return total_fit / (double) this->population_.size();
    }

    template <class IndividualType>
    double Population<IndividualType>::getPairwiseDiversity() const {
        using gene_type_t = typename std::tuple_element<0, typename IndividualType::gene_types>::type;

        double diversity = 0;

        std::vector<std::pair<gene_type_t, gene_type_t>> bounds = this->population_[0].getBounds();
        double normalization_factor = 0;
        for (size_t index = 0; index < bounds.size(); ++index)
            normalization_factor += pow(bounds[index].second - bounds[index].first, 2);
        normalization_factor = sqrt(normalization_factor);

        std::vector<double> ind_dist(this->population_.size());

#pragma omp parallel for schedule(dynamic)
        for (size_t index = 0; index < this->population_.size(); ++index) {
            for (size_t index2 = index + 1; index2 < this->population_.size(); ++index2) {
                std::vector<gene_type_t> chromosome_1 = this->population_[index].getChromosome();
                std::vector<gene_type_t> chromosome_2 = this->population_[index2].getChromosome();
                double pair_distance = 0;
                for (size_t gene_index = 0; gene_index < bounds.size(); ++gene_index) {
                    pair_distance += pow(chromosome_1[gene_index] - chromosome_2[gene_index], 2);
                }
                ind_dist[index] += sqrt(pair_distance);
            }
        }
        diversity += std::accumulate(ind_dist.begin(), ind_dist.end(), 0.0);
        diversity = (2 * diversity) / (double) (this->population_.size() * (this->population_.size() - 1));
        diversity /= normalization_factor;

        return diversity;
    }
}

#endif