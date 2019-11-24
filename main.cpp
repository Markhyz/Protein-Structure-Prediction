#include "evolutionary_algorithms/core/include/individual.hpp"
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

using namespace std;

class NullFitness : public EvoAlg::AbstractFitnessFunction {
  public:
    virtual fitness_t operator()(EvoAlg::AbstractGenotype const& genotype) const {
        return {-5.3};
    }

    virtual size_t getSize() const {
        return size;
    }

    virtual vector<int8_t> const& getSign() const {
        return sign;
    }

  private:
    size_t size = 1;
    vector<int8_t> sign{1};
};

int main() {
    shared_ptr<NullFitness> fit = make_shared<NullFitness>();
    EvoAlg::Individual<int> ind(fit, {3, 6, 2, 1, 7, 9});
    ind.setChromosome<0>(3, -10);
    for (int x : ind.getChromosome<0>()) {
        cout << x << " ";
    }
    cout << endl;
    ind.evaluateFitness();
    double x = ind.getFitnessValue()[0];
    cout << x << endl;
    ind.setChromosome<0>(vector<int>{2, 3, 5});
    for (int x : ind.getChromosome<0>()) {
        cout << x << " ";
    }
    return 0;
}