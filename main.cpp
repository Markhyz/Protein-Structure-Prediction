#include "evolutionary_algorithms/core/include/individual.hpp"
#include <iostream>
#include <queue>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

using namespace std;

int main() {
    vector<size_t> sizes = {10, 2, 31};
    vector<int> v1 = {1, 2, 5, 41};
    vector<double> v2 = {2.3, 1.002, 131.03, 0.4};
    vector<string> v3 = {"esdad", "sfddsf", "da"};
    vector<char> v4 = {'s', 'x'};
    EvoAlg::Genotype<int, double, string> genotype(sizes);
    EvoAlg::Genotype<int, double, string> genotype2(v1, v2, v3);
    for (auto &x : genotype2.getChromosome<0>()) {
        cout << x << " ";
    }
    cout << endl;
    for (auto &x : genotype2.getChromosome<1>()) {
        cout << x << " ";
    }
    cout << endl;
    for (auto &x : genotype2.getChromosome<2>()) {
        cout << x << " ";
    }
    cout << endl;
    genotype.setChromosome<2>(1, "sss");
    vector<double> vd{2.3, 0.432, 11, 234.1121};
    genotype.setChromosome<1>(vd);
    return 0;
}