#include "../../include/evo_alg/core/fitness.hpp"
#include "../../include/evo_alg/commons/utils.hpp"

namespace evo_alg {
    namespace fitness {
        std::vector<double> linearScale(std::vector<double> const& f, double const c) {
            std::vector<double> sf(f.size());
            double f_min = 1e9, f_avg = 0, f_max = -1e9, f_tot = 0;
            for (size_t i = 0; i < f.size(); ++i) {
                f_min = std::min(f_min, f[i]);
                f_max = std::max(f_max, f[i]);
                f_tot += f[i];
            }
            f_avg = f_tot / (double) f.size();
            double alfa, beta;
            if (f_min > (c * f_avg - f_max) / (c - 1)) {
                if (utils::numericEqual(f_max, f_avg)) {
                    sf = f;
                    return sf;
                }
                alfa = (f_avg * (c - 1)) / (f_max - f_avg);
                beta = (f_avg * (f_max - c * f_avg)) / (f_max - f_avg);
            } else {
                if (utils::numericEqual(f_avg, f_min)) {
                    sf = f;
                    return sf;
                }
                alfa = f_avg / (f_avg - f_min);
                beta = (-f_min * f_avg) / (f_avg - f_min);
            }

            for (size_t i = 0; i < f.size(); ++i) {
                sf[i] = f[i] * alfa + beta;
            }

            return sf;
        }

        std::vector<double> linearNormalization(std::vector<double> const& fitness_values, double const min_value,
                                                double const max_value) {
            if (utils::numericEqual(min_value, max_value))
                return fitness_values;

            std::vector<double> result(fitness_values);
            for (double& fit : result) {
                fit = (fit - min_value) / (max_value - min_value);
            }

            return result;
        }
    }
}