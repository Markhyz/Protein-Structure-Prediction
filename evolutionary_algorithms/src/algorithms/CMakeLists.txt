
add_library(ga STATIC ga.cpp)
target_link_libraries(ga PUBLIC evo_alg)
target_compile_features(ga PUBLIC cxx_std_17)
target_compile_options(ga PRIVATE ${CXX_FLAGS})

add_library(algorithms INTERFACE)
target_link_libraries(algorithms INTERFACE ga)
target_compile_features(algorithms INTERFACE cxx_std_17)

