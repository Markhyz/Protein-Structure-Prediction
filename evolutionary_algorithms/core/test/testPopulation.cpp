#include "../include/fitness.hpp"
#include "../include/individual.hpp"
#include "../include/population.hpp"

#include <memory>
#include <numeric>

#include "gtest/gtest.h"

using namespace EvoAlg;

class SimpleFitness : public AbstractFitnessFunction<int> {
  public:
    virtual fitness_t operator()(Genotype<int> const& genotype) const {
        std::vector<int> v = genotype.getChromosome();

        return {std::accumulate(v.begin(), v.end(), 0.0)};
    }

    virtual size_t getDimension() const {
        return 1;
    }

    virtual std::vector<bool> const& getDirection() const {
        return direction;
    }

    std::vector<bool> direction{minimize};
};

class PopulationTest : public ::testing::Test {
  protected:
    SimpleFitness::shared_ptr f{std::make_shared<SimpleFitness>()};
    Individual<int> individual_1{f, {10, 20, 123, 40}};
    Individual<int> individual_2{f, {0, -5, 13}};
    Individual<int> individual_3{f, {1111, 10000, 200, 600, 50}};
    Population<Individual<int>> population{{individual_1, individual_2, individual_3}};
};

TEST_F(PopulationTest, Initialize) {
    EXPECT_EQ(individual_1.getChromosome(), population.getIndividual(0).getChromosome());
    EXPECT_EQ(individual_2.getChromosome(), population.getIndividual(1).getChromosome());
    EXPECT_EQ(individual_3.getChromosome(), population.getIndividual(2).getChromosome());

    EXPECT_EQ(individual_1.getFitnessFunction(), population.getIndividual(0).getFitnessFunction());
    EXPECT_EQ(individual_2.getFitnessFunction(), population.getIndividual(1).getFitnessFunction());
    EXPECT_EQ(individual_3.getFitnessFunction(), population.getIndividual(2).getFitnessFunction());
}

TEST_F(PopulationTest, EvaluateFitness) {
    population.evaluateFitness();
    int ind1_f = (int) population.getIndividual(0).getFitnessValue()[0];
    int ind2_f = (int) population.getIndividual(1).getFitnessValue()[0];
    int ind3_f = (int) population.getIndividual(2).getFitnessValue()[0];

    EXPECT_EQ(193, ind1_f);
    EXPECT_EQ(8, ind2_f);
    EXPECT_EQ(11961, ind3_f);

    Individual<int> new_ind{f, {241, -500}};
    new_ind.evaluateFitness();

    population.setIndividual(1, new_ind);

    ind2_f = (int) population.getIndividual(1).getFitnessValue()[0];

    EXPECT_EQ(-259, ind2_f);
}

TEST_F(PopulationTest, addNewIndividual) {
    size_t last_idx = population.getSize() - 1;

    EXPECT_EQ(individual_3.getChromosome(), population.getIndividual(last_idx).getChromosome());
    EXPECT_EQ(individual_3.getFitnessFunction(), population.getIndividual(last_idx).getFitnessFunction());

    Individual<int> new_ind{f, {100, 100, 1000}};
    population.appendIndividual(new_ind);
    last_idx = population.getSize() - 1;

    EXPECT_EQ(new_ind.getChromosome(), population.getIndividual(last_idx).getChromosome());
    EXPECT_EQ(new_ind.getFitnessFunction(), population.getIndividual(last_idx).getFitnessFunction());
}