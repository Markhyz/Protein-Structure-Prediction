#include "../../include/evo_alg/commons/utils.hpp"

namespace evo_alg {
    namespace utils {
        std::random_device rd;
        std::mt19937_64 rng(rd());

        double uniformProbGen() {
            std::uniform_real_distribution<double> uniform_dist(0, 1);
            return uniform_dist(rng);
        }
    }
}