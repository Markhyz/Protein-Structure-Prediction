#include "../include/fitness.hpp"
#include "../include/individual.hpp"

#include <memory>
#include <numeric>

#include "gtest/gtest.h"

using namespace EvoAlg;

class SimpleFitness : public AbstractFitnessFunction {
  public:
    virtual fitness_t operator()(AbstractGenotype const& genotype) const {
        std::vector<double> chromosome = dynamic_cast<Genotype<double> const&>(genotype).getChromosome();
        double total = std::accumulate(chromosome.begin(), chromosome.end(), 0.0);

        return {total};
    }

    virtual size_t getDimension() const {
        return dimension;
    }

    virtual std::vector<bool> const& getDirection() const {
        return sign;
    }

  private:
    size_t dimension = 1;
    std::vector<bool> sign{minimize};
};

class IndividualTest : public ::testing::Test {
  protected:
    IndividualTest() {
        f = std::make_shared<SimpleFitness>();
        chromosome = {2.7, 5.3, -15.2, 11, 57, 1.062};
        individual = Individual<double>(f, chromosome);
    }

    std::shared_ptr<SimpleFitness> f;
    std::vector<double> chromosome;
    Individual<double> individual;
};

TEST_F(IndividualTest, Initialize) {
    std::vector<double> ind_chromosome = individual.getChromosome();

    EXPECT_EQ(chromosome.size(), ind_chromosome.size());

    for (size_t i = 0; i < chromosome.size(); ++i)
        EXPECT_DOUBLE_EQ(chromosome[i], ind_chromosome[i]);
}

TEST_F(IndividualTest, EvaluateFitness) {
    individual.evaluateFitness();
    double total = individual.getFitnessValue()[0];

    EXPECT_DOUBLE_EQ(61.862, total);
}

class CompositeFitness : public AbstractFitnessFunction {
  public:
    virtual fitness_t operator()(AbstractGenotype const& genotype) const {
        std::vector<double> chromosome_1 =
            dynamic_cast<Genotype<double, bool, char> const&>(genotype).getChromosome<0>();
        std::vector<bool> chromosome_2 = dynamic_cast<Genotype<double, bool, char> const&>(genotype).getChromosome<1>();
        std::vector<char> chromosome_3 = dynamic_cast<Genotype<double, bool, char> const&>(genotype).getChromosome<2>();

        double total = std::accumulate(chromosome_1.begin(), chromosome_1.end(), 0.0);
        int x =
            std::accumulate(chromosome_2.begin(), chromosome_2.end(), 0, [](int tot, bool cur) { return tot + cur; });
        int y = std::accumulate(chromosome_3.begin(), chromosome_3.end(), 0,
                                [](int tot, char cur) { return tot + (cur == 'a' || cur == 't'); });

        return {total, (double)x * y};
    }

    virtual size_t getDimension() const {
        return dimension;
    }

    virtual std::vector<bool> const& getDirection() const {
        return sign;
    }

  private:
    size_t dimension = 1;
    std::vector<bool> sign{minimize};
};

class CompositeIndividualTest : public ::testing::Test {
  protected:
    CompositeIndividualTest() {
        f = std::make_shared<CompositeFitness>();
        chromosome_1 = {5.43, -11, 103.4};
        chromosome_2 = {0, 0, 1, 1, 0, 1};
        chromosome_3 = {'a', 't', 't', 'g', 'c', 'a', 't', 'a'};
        individual = Individual<double, bool, char>(f, chromosome_1, chromosome_2, chromosome_3);
    }

    std::shared_ptr<CompositeFitness> f;
    std::vector<double> chromosome_1;
    std::vector<bool> chromosome_2;
    std::vector<char> chromosome_3;
    Individual<double, bool, char> individual;
};

TEST_F(CompositeIndividualTest, Initialize) {
    std::vector<double> ind_chromosome_1 = individual.getChromosome<0>();
    std::vector<bool> ind_chromosome_2 = individual.getChromosome<1>();
    std::vector<char> ind_chromosome_3 = individual.getChromosome<2>();

    ASSERT_EQ(chromosome_1.size(), ind_chromosome_1.size());
    for (size_t i = 0; i < chromosome_1.size(); ++i)
        EXPECT_DOUBLE_EQ(chromosome_1[i], ind_chromosome_1[i]);

    ASSERT_EQ(chromosome_2.size(), ind_chromosome_2.size());
    for (size_t i = 0; i < chromosome_2.size(); ++i)
        EXPECT_EQ(chromosome_2[i], ind_chromosome_2[i]);

    ASSERT_EQ(chromosome_3.size(), ind_chromosome_3.size());
    for (size_t i = 0; i < chromosome_3.size(); ++i)
        EXPECT_EQ(chromosome_3[i], ind_chromosome_3[i]);
}

TEST_F(CompositeIndividualTest, EvaluateFitness) {
    individual.evaluateFitness();
    std::vector<double> total = individual.getFitnessValue();

    EXPECT_DOUBLE_EQ(97.83, total[0]);
    EXPECT_EQ(18, total[1]);
}

TEST_F(CompositeIndividualTest, ChangeChromosome) {
    std::vector<double> new_chromosome_1 = {10, 20, 30};
    std::vector<bool> new_chromosome_2 = {1, 1, 1, 1, 1};

    individual.setChromosome<1>(new_chromosome_2);

    std::vector<bool> ind_chromosome_2 = individual.getChromosome<1>();

    ASSERT_EQ(new_chromosome_2.size(), ind_chromosome_2.size());
    for (size_t i = 0; i < chromosome_2.size(); ++i)
        EXPECT_EQ(new_chromosome_2[i], ind_chromosome_2[i]);

    individual.setChromosome<0>(new_chromosome_1);

    std::vector<double> ind_chromosome_1 = individual.getChromosome<0>();

    ASSERT_EQ(new_chromosome_1.size(), ind_chromosome_1.size());
    for (size_t i = 0; i < chromosome_1.size(); ++i)
        EXPECT_DOUBLE_EQ(new_chromosome_1[i], ind_chromosome_1[i]);

    individual.evaluateFitness();
    std::vector<double> total = individual.getFitnessValue();

    EXPECT_DOUBLE_EQ(60, total[0]);
    EXPECT_EQ(30, total[1]);
}