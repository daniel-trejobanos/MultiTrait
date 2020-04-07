set(BINARY ${CMAKE_PROJECT_NAME})


add_library(covsampler
  GammaProp.cpp
  GammaProp.h
  NormalProp.cpp
  NormalProp.h
  SampleVariance.cpp
  TargetDist.h
  VarianceTarget.cpp
  VarianceTarget.h
  gammaNewtonian.h
  gaussianNewtonian.h
  groupVarianceTarget.cpp
  groupVarianceTarget.h
  newtonian.cpp
  newtonian.h
  proposal.h
  sampleGroupVariances.cpp
)

add_executable(${BINARY}_run main.cpp)


target_link_libraries(${BINARY}_run covsampler)


set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS ${CMMAKE_PROJECT_SOURCE_DIR} ${CMAKE_PROJECT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} CACHE INTERNAL "${CMAKE_PROJECT_NAME}: Include Directories"
    FORCE)
