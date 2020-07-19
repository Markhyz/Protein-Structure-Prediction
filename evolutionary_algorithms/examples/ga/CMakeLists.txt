add_executable(integer_sum integer_sum.cpp)
target_link_libraries(integer_sum PRIVATE evo_alg)
target_compile_features(integer_sum PRIVATE cxx_std_17)
target_compile_options(integer_sum PRIVATE ${CXX_FLAGS})