add_library(sbx STATIC sbx.cpp)
target_link_libraries(sbx PUBLIC evo_alg)
target_compile_features(sbx PUBLIC cxx_std_17)
target_compile_options(sbx PRIVATE ${CXX_FLAGS})

add_library(point STATIC point.cpp)
target_link_libraries(point PUBLIC evo_alg)
target_compile_features(point PUBLIC cxx_std_17)
target_compile_options(point PRIVATE ${CXX_FLAGS})

add_library(recombination INTERFACE)
target_link_libraries(recombination INTERFACE sbx point)
target_compile_features(recombination INTERFACE cxx_std_17)
