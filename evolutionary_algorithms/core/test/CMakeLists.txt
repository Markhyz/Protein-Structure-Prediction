set(CXX_FLAGS -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wunused -Wconversion -Wsign-conversion -g)

add_executable(test_individual ./testIndividual.cpp)
target_include_directories(test_individual PRIVATE ../include)
target_link_libraries(test_individual PRIVATE individual gtest_main)
target_compile_features(test_individual PRIVATE cxx_std_17)
target_compile_options(test_individual PRIVATE ${CXX_FLAGS})

add_executable(test_population ./testPopulation.cpp)
target_include_directories(test_population PRIVATE ../include)
target_link_libraries(test_population PRIVATE population gtest_main)
target_compile_features(test_population PRIVATE cxx_std_17)
target_compile_options(test_population PRIVATE ${CXX_FLAGS})