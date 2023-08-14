cmake_minimum_required(VERSION 3.17)

project(test LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(evolutionary_algorithms)
add_subdirectory(external/rosetta)

add_executable(mobrkga_pred main.cpp)
target_link_libraries(mobrkga_pred PRIVATE rosetta evo_alg)
target_compile_features(mobrkga_pred PUBLIC cxx_std_17)
target_compile_options(mobrkga_pred PRIVATE -g)

add_executable(score_calc scripts/score_calc.cpp)
target_link_libraries(score_calc PRIVATE rosetta)
target_compile_features(score_calc PUBLIC cxx_std_17)

add_executable(population_score scripts/population_score.cpp)
target_link_libraries(population_score PRIVATE rosetta)
target_compile_features(population_score PUBLIC cxx_std_17)

if (OMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(population_score PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

if (MPI)
    find_package(MPI REQUIRED)
    message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

    target_link_libraries(mobrkga_pred PRIVATE MPI::MPI_CXX)
    target_compile_options(mobrkga_pred PRIVATE -DEVO_ALG_MPI)
endif()
