#ifndef GUARD_H_FITNESS
#define GUARD_H_FITNESS

#include "../../../commons/include/macros.hpp"
#include "../../../commons/include/types.hpp"
#include "genotype.hpp"
#include <vector>

namespace EvoAlg {
    class AbstractFitness {
      public:
        POINTER_ALIAS(AbstractFitness)

        using fitness_t = std::vector<double>;

        virtual fitness_t operator()(AbstractGenotype const& genotype) const = 0;

        size_t getSize() const;
        int8_t getSign(size_t index) const;

        virtual ~AbstractFitness() = 0;

      private:
        size_t size_;
        std::vector<int8_t> sign_;
    };
}

#endif