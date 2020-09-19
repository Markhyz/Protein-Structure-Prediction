cmake_minimum_required(VERSION 3.17)

project(test LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(evolutionary_algorithms)
add_subdirectory(external/rosetta)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE rosetta evo_alg)
target_compile_features(main PUBLIC cxx_std_17)

if (MPI)
    find_package(MPI REQUIRED)
    message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
    target_link_libraries(main PRIVATE MPI::MPI_CXX)
    target_compile_options(main PRIVATE -DEVO_ALG_MPI)
endif()