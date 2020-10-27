#include "../../include/evo_alg/commons/utils.hpp"

namespace evo_alg {
    namespace utils {
        std::random_device rd;
        std::mt19937_64 rng(rd());

        int mpi_size = 1, mpi_rank = 1;

        double uniformProbGen() {
            std::uniform_real_distribution<double> uniform_dist(0, 1);
            return uniform_dist(rng);
        }

        void initMPI(int const size, int const rank) {
            mpi_size = size;
            mpi_rank = rank;
        }
    }
}