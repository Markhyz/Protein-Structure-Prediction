cmake_minimum_required(VERSION 3.17)

project(EvolutionaryAlgorithms VERSION 1.0.0
                               DESCRIPTION "Evolutionary Algorithms library"
                               LANGUAGES CXX)

add_library(evo_alg INTERFACE)
target_include_directories(evo_alg INTERFACE include)
target_compile_features(evo_alg INTERFACE cxx_std_17)

add_subdirectory(src)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(tests)
endif()