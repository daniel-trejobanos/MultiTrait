cmake_minimum_required(VERSION 3.10)
project(mtrait)
set(CMAKE_PROJECT_NAME mtrait)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
if (NOT TARGET Eigen3::Eigen)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/Eigen/)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen/)
    SET(GCC_EIGEN_COMPILE_FLAGS "-Wno-ignored-attributes -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DEIGEN_UNROLLING_LIMIT=1000")
  else()
    message(FATAL_ERROR "Cound not find Eigen submodule!")
  endif()
endif()


Find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()
include(CheckCXXCompilerFlag)

# Check for -march=native support in the compiler
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()


find_package(Boost REQUIRED COMPONENTS locale)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    message("-- Boost found: ${Boost_VERSION_STRING}")
else()
    message(FATAL_ERROR "Boost not found!")
endif()


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_EIGEN_COMPILE_FLAGS}")


include_directories(src)

add_subdirectory(src)


if (CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR
    CMAKE_BUILD_TYPE STREQUAL "RelWithDebugInfo" )
  add_definitions("-DDEBUG")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  
elseif ( CMAKE_BUILD_TYPE STREQUAL "RELEASE" )
    add_definitions("-O3")
endif()




if(ENABLE_UNIT_TESTS)
  if(NOT TARGET gtest)  
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
      )
    execute_process(COMMAND "${CMAKE_COMMAND}" --build .
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
      )

    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This adds the following targets:
    # gtest, gtest_main, gmock and gmock_main
    add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
      "${CMAKE_BINARY_DIR}/googletest-build"
      )

    # The gtest/gmock targets carry header search path dependencies
    # automatically when using CMake 2.8.11 or later. Otherwise we
    # have to add them here ourselves.
    if(CMAKE_VERSION VERSION_LESS 2.8.11)
      include_directories("${gtest_SOURCE_DIR}/include"
	"${gmock_SOURCE_DIR}/include"
	)
    endif()

    # Now simply link your own targets against gtest, gmock,
    # etc. as appropriate

    add_subdirectory(tst)
  endif()
endif()
