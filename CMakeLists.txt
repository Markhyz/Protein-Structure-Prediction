cmake_minimum_required(VERSION 3.17)

project(test LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(evolutionary_algorithms)
add_subdirectory(external/rosetta)

add_executable(ga_main ga_main.cpp)
target_link_libraries(ga_main PRIVATE rosetta evo_alg)
target_compile_features(ga_main PUBLIC cxx_std_17)

add_executable(brkga_main brkga_main.cpp)
target_link_libraries(brkga_main PRIVATE rosetta evo_alg)
target_compile_features(brkga_main PUBLIC cxx_std_17)

add_executable(mobrkga_main mobrkga_main.cpp)
target_link_libraries(mobrkga_main PRIVATE rosetta evo_alg)
target_compile_features(mobrkga_main PUBLIC cxx_std_17)

add_executable(score_calc score_calc.cpp)
target_link_libraries(score_calc PRIVATE rosetta)
target_compile_features(score_calc PUBLIC cxx_std_17)

if (MPI)
    find_package(MPI REQUIRED)
    message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
    target_link_libraries(ga_main PRIVATE MPI::MPI_CXX)
    target_compile_options(ga_main PRIVATE -DEVO_ALG_MPI)

    target_link_libraries(brkga_main PRIVATE MPI::MPI_CXX)
    target_compile_options(brkga_main PRIVATE -DEVO_ALG_MPI)
endif()