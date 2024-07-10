cmake_minimum_required(VERSION 3.14...3.30)

project(factorization
        LANGUAGES CXX
        DESCRIPTION "Sandbox to play with factorization")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python 3.12 COMPONENTS Interpreter Development)

execute_process(COMMAND python -c "import numpy; print(numpy.get_include())"
                RESULT_VARIABLE result
                OUTPUT_VARIABLE NUMPY_INCLUDE)

if (result EQUAL 0)
  string(STRIP ${NUMPY_INCLUDE} NUMPY_INCLUDE)
  message("Numpy include directory: ${NUMPY_INCLUDE}")
else()
  message(FATAL_ERROR "Could not find numpy include directory")
endif()

# Begin of gitache configuration.
set(GITACHE_PACKAGES libcwd_r)
include(cwm4/cmake/StableGitache)
# End of gitache configuration.

# This project uses aicxx modules.
include(cwm4/cmake/AICxxProject)

#Exit if someone tries to contaminates the source directory with an in-source build.
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "Please do out-of-source builds.\nCleanup: \"rm -rf CMake CMakeCache.txt CMakeFiles/\"")
endif()

include(AICxxSubmodules)

add_subdirectory(src)
