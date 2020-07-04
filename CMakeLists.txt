cmake_minimum_required(VERSION 3.17)

project(test LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(evolutionary_algorithms)

add_executable(test main.cpp)
target_link_libraries(test PUBLIC individual)
target_compile_features(test PRIVATE cxx_std_17)