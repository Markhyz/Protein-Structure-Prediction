cmake_minimum_required(VERSION 3.17)

project(test LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(evolutionary_algorithms)
add_subdirectory(external/rosetta)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE rosetta evo_alg)
target_compile_features(main PUBLIC cxx_std_17)