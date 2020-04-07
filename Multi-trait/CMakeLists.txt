cmake_minimum_required(VERSION 3.10)
project(mtrait)

set(CMAKE_CXX_STANDARD 14)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/Eigen/)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen/)
  SET(GCC_EIGEN_COMPILE_FLAGS "-Wno-ignored-attributes -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DEIGEN_UNROLLING_LIMIT=1000")
else()
    message(FATAL_ERROR "Cound not find Eigen submodule!")
endif()


find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()
include(CheckCXXCompilerFlag)

# Check for -march=native support in the compiler
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()



SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_EIGEN_COMPILE_FLAGS}")


include_directories(src)

add_subdirectory(src)
add_subdirectory(lib/googletest)

if(ENABLE_UNIT_TESTS)
    add_subdirectory(tst)
endif()
find_package(Boost REQUIRED COMPONENTS locale)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    message("-- Boost found: ${Boost_VERSION_STRING}")
else()
    message(FATAL_ERROR "Boost not found!")
endif()
