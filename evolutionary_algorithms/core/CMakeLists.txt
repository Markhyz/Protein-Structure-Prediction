set(CXX_FLAGS -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wunused -Wconversion -Wsign-conversion)

add_library(fitness INTERFACE)
target_include_directories(fitness INTERFACE ./include)
target_link_libraries(fitness INTERFACE commons)
target_compile_features(fitness INTERFACE cxx_std_17)

add_library(phenotype STATIC ./src/phenotype.cpp)
target_include_directories(phenotype PUBLIC ./include)
target_link_libraries(phenotype PUBLIC fitness commons)
target_compile_features(phenotype PUBLIC cxx_std_17)
target_compile_options(phenotype PRIVATE ${CXX_FLAGS})

add_library(genotype STATIC ./src/genotype.cpp)
target_include_directories(genotype PUBLIC ./include)
target_link_libraries(genotype PUBLIC commons)
target_compile_features(genotype PUBLIC cxx_std_17)
target_compile_options(genotype PRIVATE ${CXX_FLAGS})

add_library(individual STATIC ./src/individual.cpp)
target_include_directories(individual PUBLIC ./include)
target_link_libraries(individual PUBLIC phenotype genotype commons)
target_compile_features(individual PUBLIC cxx_std_17)
target_compile_options(individual PRIVATE ${CXX_FLAGS})


