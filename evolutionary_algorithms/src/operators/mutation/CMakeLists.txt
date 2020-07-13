add_library(bitFlip STATIC bitFlip.cpp)
target_link_libraries(bitFlip PUBLIC evo_alg)
target_compile_features(bitFlip PUBLIC cxx_std_17)
target_compile_options(bitFlip PRIVATE ${CXX_FLAGS})

add_library(polynomial STATIC polynomial.cpp)
target_link_libraries(polynomial PUBLIC evo_alg)
target_compile_features(polynomial PUBLIC cxx_std_17)
target_compile_options(polynomial PRIVATE ${CXX_FLAGS})

add_library(mutation INTERFACE)
target_link_libraries(mutation INTERFACE bitFlip polynomial)
target_compile_features(mutation INTERFACE cxx_std_17)
