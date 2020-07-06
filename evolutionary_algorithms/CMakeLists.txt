cmake_minimum_required(VERSION 3.17)

project(EvolutionaryAlgorithms VERSION 1.0.0
                               DESCRIPTION "Evolutionary Algorithms library"
                               LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)

    include(FetchContent)

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        release-1.10.0
    )

    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(commons)
add_subdirectory(core)
