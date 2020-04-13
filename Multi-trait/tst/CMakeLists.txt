set(BINARY ${CMAKE_PROJECT_NAME}_tst)

add_executable(${BINARY}
  linModTst.cpp)


set_property(TARGET ${BINARY} PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ${BINARY} PROPERTY CXX_STANDARD 17)


target_link_libraries(${BINARY} ${RCPPLIBS_l})
target_link_libraries(${BINARY} ${RINSIDELIBS_l})
target_link_libraries(${BINARY} ${RLDFLAGS_l})
target_link_libraries(${BINARY}  ${CMAKE_PROJECT_NAME}_lib gtest)


add_test(NAME ${BINARY} COMMAND ${BINARY})
