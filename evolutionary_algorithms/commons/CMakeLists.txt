

add_library(commons INTERFACE)
target_include_directories(commons INTERFACE ./include)
target_compile_features(commons INTERFACE cxx_std_17)
