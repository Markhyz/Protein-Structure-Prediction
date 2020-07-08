set(CXX_FLAGS -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wunused -Wconversion -Wsign-conversion)

add_library(fitness INTERFACE)
target_include_directories(fitness INTERFACE ./include)
target_link_libraries(fitness INTERFACE commons)
target_compile_features(fitness INTERFACE cxx_std_17)

add_library(phenotype INTERFACE)
target_include_directories(phenotype INTERFACE ./include)
target_link_libraries(phenotype INTERFACE fitness commons)
target_compile_features(phenotype INTERFACE cxx_std_17)

add_library(genotype INTERFACE)
target_include_directories(genotype INTERFACE ./include)
target_link_libraries(genotype INTERFACE commons)
target_compile_features(genotype INTERFACE cxx_std_17)

add_library(individual INTERFACE)
target_include_directories(individual INTERFACE ./include)
target_link_libraries(individual INTERFACE phenotype genotype commons)
target_compile_features(individual INTERFACE cxx_std_17)

add_library(population INTERFACE)
target_include_directories(population INTERFACE ./include)
target_link_libraries(population INTERFACE individual commons)
target_compile_features(population INTERFACE cxx_std_17)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(test)
endif()

