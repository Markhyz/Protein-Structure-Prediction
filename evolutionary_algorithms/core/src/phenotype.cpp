#include "../include/phenotype.hpp"
#include <iostream>

using namespace EvoAlg;

Phenotype::Phenotype(){};

Phenotype::Phenotype(typename AbstractFitnessFunction::const_shared_ptr const& fitness) : fitness_(fitness){};

void Phenotype::evaluateFitness(AbstractGenotype const& genotype) {
    fitness_value_ = (*fitness_)(genotype);
}

typename AbstractFitnessFunction::const_shared_ptr Phenotype::getFitnessFunction() const {
    return fitness_;
}

typename AbstractFitnessFunction::fitness_t const& Phenotype::getFitnessValue() const {
    if (!fitness_value_.has_value()) {
        throw UndefinedFitnessException();
    }

    return fitness_value_.value();
}