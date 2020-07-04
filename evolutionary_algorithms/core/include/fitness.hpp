#ifndef GUARD_H_FITNESS
#define GUARD_H_FITNESS

#include "../../commons/include/macros.hpp"
#include "../../commons/include/types.hpp"
#include "genotype.hpp"
#include <vector>

namespace EvoAlg {
    class AbstractFitnessFunction {
      public:
        POINTER_ALIAS(AbstractFitnessFunction)

        using fitness_t = std::vector<double>;

        virtual ~AbstractFitnessFunction() = default;

        virtual std::vector<bool> const& getDirection() const = 0;
        virtual size_t getDimension() const = 0;

        virtual fitness_t operator()(AbstractGenotype const& genotype) const = 0;

        static constexpr bool minimize = 0;
        static constexpr bool maximize = 1;
    };
}

#endif