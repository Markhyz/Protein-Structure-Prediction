#ifndef GUARD_H_FITNESS
#define GUARD_H_FITNESS

#include "../../../commons/include/macros.hpp"
#include "../../../commons/include/types.hpp"
#include "genotype.hpp"
#include <vector>

namespace EvoAlg {
    class AbstractFitnessFunction {
      public:
        POINTER_ALIAS(AbstractFitnessFunction)

        using fitness_t = std::vector<double>;

        virtual std::vector<int8_t> const& getSign() const = 0;
        virtual size_t getSize() const = 0;

        virtual fitness_t operator()(AbstractGenotype const& genotype) const = 0;

        virtual ~AbstractFitnessFunction() = default;
    };
}

#endif