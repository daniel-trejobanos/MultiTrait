set(BINARY ${CMAKE_PROJECT_NAME})


add_library(${CMAKE_PROJECT_NAME}_lib STATIC
  GammaProp.cpp
  GammaProp.h
  NormalProp.cpp
  NormalProp.h
  SampleVariance.cpp
  TargetDist.h
  VarianceTarget.cpp
  VarianceTarget.h
  groupVarianceTarget.cpp
  groupVarianceTarget.h
  newtonian.cpp
  newtonian.h
  proposal.h
  sampleGroupVariances.cpp
  linearModel.h
  linearModel.cpp
  dist.h
  dist.cpp
  epsilonVarianceTarget.cpp
  epsilonVarianceTarget.h
  sampleEpsVariance.cpp
  sampleEpsVariance.h
  InvGammaProp.h
  InvGammaProp.cpp
)

add_executable(${BINARY}_run main.cpp)


target_link_libraries(${BINARY}_run ${CMAKE_PROJECT_NAME}_lib) 

Set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS ${CMMAKE_PROJECT_SOURCE_DIR} ${CMAKE_PROJECT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} CACHE INTERNAL "${CMAKE_PROJECT_NAME}: Include Directories"
    FORCE)
