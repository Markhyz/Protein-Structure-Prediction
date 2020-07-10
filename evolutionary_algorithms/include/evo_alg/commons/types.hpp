#ifndef GUARD_H_EVO_ALG_TYPES
#define GUARD_H_EVO_ALG_TYPES

#include <cstddef>
#include <tuple>

template <std::size_t N, typename... Ts>
using NthType = typename std::tuple_element<N, std::tuple<Ts...>>::type;

#endif