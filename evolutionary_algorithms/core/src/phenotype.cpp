#include "../include/phenotype.hpp"
#include <stdexcept>

using namespace EvoAlg;

Phenotype::Phenotype(AbstractFitness::shared_ptr const& fitness) : fitness_(fitness){};

double Phenotype::getFitnessValue(size_t index) const {
    if (!fitness_value_.has_value()) {
        throw Phenotype::UndefinedFitnessException();
    }

    std::vector<double> const& value = fitness_value_.value();
    if (index < 0 || index >= value.size()) {
        throw std::out_of_range("out of range access");
    }

    return value[index];
}

void Phenotype::evaluateFitness(AbstractGenotype const& genotype) {
    fitness_value_ = (*fitness_)(genotype);
}

const char* Phenotype::UndefinedFitnessException::what() const throw() {
    return "undefined fitness value";
}