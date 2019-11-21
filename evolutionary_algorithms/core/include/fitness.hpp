#ifndef GUARD_H_FITNESS
#define GUARD_H_FITNESS

#include "../../../commons/include/macros.hpp"
#include "../../../commons/include/types.hpp"
#include "genotype.hpp"
#include <vector>

namespace EvoAlg {
    template <size_t FitnessSize>
    class AbstractFitnessFunction {
      public:
        POINTER_ALIAS(AbstractFitnessFunction<FitnessSize>)

        using fitness_t = std::conditional_t<FitnessSize == 1, double, std::vector<double>>;

        virtual int8_t getSign(size_t index) const = 0;

        virtual fitness_t operator()(AbstractGenotype const& genotype) const = 0;

        virtual ~AbstractFitnessFunction() = default;

        constexpr static size_t size = FitnessSize;
    };
}

#endif