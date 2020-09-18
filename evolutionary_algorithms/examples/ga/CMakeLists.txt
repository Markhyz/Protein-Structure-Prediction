add_executable(integer_sum integer_sum.cpp)
target_link_libraries(integer_sum PRIVATE evo_alg)
target_compile_features(integer_sum PRIVATE cxx_std_17)
target_compile_options(integer_sum PRIVATE ${CXX_FLAGS})

add_executable(sphere sphere.cpp)
target_link_libraries(sphere PRIVATE evo_alg)
target_compile_features(sphere PRIVATE cxx_std_17)
target_compile_options(sphere PRIVATE ${CXX_FLAGS})

add_executable(ackley ackley.cpp)
target_link_libraries(ackley PRIVATE evo_alg)
target_compile_features(ackley PRIVATE cxx_std_17)
target_compile_options(ackley PRIVATE ${CXX_FLAGS})

add_executable(step step.cpp)
target_link_libraries(step PRIVATE evo_alg)
target_compile_features(step PRIVATE cxx_std_17)
target_compile_options(step PRIVATE ${CXX_FLAGS})

add_executable(gp1 gp1.cpp)
target_link_libraries(gp1 PRIVATE evo_alg)
target_compile_features(gp1 PRIVATE cxx_std_17)
target_compile_options(gp1 PRIVATE ${CXX_FLAGS})

add_executable(schaffer schaffer.cpp)
target_link_libraries(schaffer PRIVATE evo_alg)
target_compile_features(schaffer PRIVATE cxx_std_17)
target_compile_options(schaffer PRIVATE ${CXX_FLAGS})

add_executable(rosenbrock rosenbrock.cpp)
target_link_libraries(rosenbrock PRIVATE evo_alg)
target_compile_features(rosenbrock PRIVATE cxx_std_17)
target_compile_options(rosenbrock PRIVATE ${CXX_FLAGS})

add_executable(schwefel schwefel.cpp)
target_link_libraries(schwefel PRIVATE evo_alg)
target_compile_features(schwefel PRIVATE cxx_std_17)
target_compile_options(schwefel PRIVATE ${CXX_FLAGS})