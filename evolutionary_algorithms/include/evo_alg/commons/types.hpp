#ifndef GUARD_H_TYPES
#define GUARD_H_TYPES

#include <cstddef>
#include <tuple>

template <std::size_t N, typename... Ts>
using NthType = typename std::tuple_element<N, std::tuple<Ts...>>::type;

#endif