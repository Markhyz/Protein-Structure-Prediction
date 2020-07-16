#ifndef GUARD_H_EVO_ALG_FITNESS
#define GUARD_H_EVO_ALG_FITNESS

#include "../commons/macros.hpp"
#include "../commons/types.hpp"

#include "genotype.hpp"

#include <vector>

namespace evo_alg {
    template <typename... ChromosomeTypes>
    class FitnessFunction {
      public:
        POINTER_ALIAS(FitnessFunction<ChromosomeTypes...>)

        using fitness_t = std::vector<double>;

        using gene_bounds_t = std::variant<std::pair<ChromosomeTypes, ChromosomeTypes>...>;
        using chromosome_bounds_t = std::vector<gene_bounds_t>;

        FitnessFunction();
        FitnessFunction(std::vector<std::pair<ChromosomeTypes, ChromosomeTypes>> const&... bounds);

        virtual ~FitnessFunction() = default;

        virtual std::vector<bool> const& getDirection() const = 0;
        virtual size_t getDimension() const = 0;

        template <size_t ChromosomeIndex = 0>
        std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                              types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
        getBounds() const;

        void setBounds(std::vector<std::pair<ChromosomeTypes, ChromosomeTypes>> const&... bounds);

        template <size_t ChromosomeIndex = 0>
        void
        setBounds(std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                                        types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const& chromosome_bounds);

        template <size_t ChromosomeIndex = 0>
        void setBounds(size_t const gene_index,
                       std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                                 types::NthType<ChromosomeIndex, ChromosomeTypes...>> const gene_bounds);

        virtual FitnessFunction* clone() const = 0;

        virtual fitness_t operator()(Genotype<ChromosomeTypes...> const& genotype) const = 0;

        static constexpr bool minimize = 0;
        static constexpr bool maximize = 1;

      private:
        std::vector<chromosome_bounds_t> bounds_;
    };

    template <typename... ChromosomeTypes>
    FitnessFunction<ChromosomeTypes...>::FitnessFunction(){};

    template <typename... ChromosomeTypes>
    FitnessFunction<ChromosomeTypes...>::FitnessFunction(
        std::vector<std::pair<ChromosomeTypes, ChromosomeTypes>> const&... bounds)
        : bounds_{
              typename FitnessFunction<ChromosomeTypes...>::chromosome_bounds_t(bounds.begin(), bounds.end())...} {};

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                          types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const
    FitnessFunction<ChromosomeTypes...>::getBounds() const {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                              types::NthType<ChromosomeIndex, ChromosomeTypes...>>>
            bounds;
        std::transform(bounds_[ChromosomeIndex].begin(), bounds_[ChromosomeIndex].end(), std::back_inserter(bounds),
                       [](gene_bounds_t const& gene_bounds) { return std::get<ChromosomeIndex>(gene_bounds); });

        return bounds;
    }

    template <typename... ChromosomeTypes>
    void FitnessFunction<ChromosomeTypes...>::setBounds(
        std::vector<std::pair<ChromosomeTypes, ChromosomeTypes>> const&... bounds) {
        bounds_ = {typename FitnessFunction<ChromosomeTypes...>::chromosome_bounds_t(bounds.begin(), bounds.end())...};
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void FitnessFunction<ChromosomeTypes...>::setBounds(
        std::vector<std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                              types::NthType<ChromosomeIndex, ChromosomeTypes...>>> const& chromosome_bounds) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));

        bounds_[ChromosomeIndex] = chromosome_bounds_t(chromosome_bounds.begin(), chromosome_bounds.end());
    }

    template <typename... ChromosomeTypes>
    template <size_t ChromosomeIndex>
    void FitnessFunction<ChromosomeTypes...>::setBounds(
        size_t const gene_index, std::pair<types::NthType<ChromosomeIndex, ChromosomeTypes...>,
                                           types::NthType<ChromosomeIndex, ChromosomeTypes...>> const gene_bounds) {
        static_assert(ChromosomeIndex < sizeof...(ChromosomeTypes));
        if (gene_index >= bounds_[ChromosomeIndex].size()) {
            throw std::out_of_range("out of range access");
        }

        bounds_[ChromosomeIndex][gene_index] = gene_bounds;
    }
}

#endif