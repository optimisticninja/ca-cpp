project(tests NONE)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})
add_compile_definitions(DEBUG)
# FIXME: Can we grab objects from other directory? Gross. Also does it fix debugging
add_executable(tests main.cc ca/ca1d_test.cc ../src/util/states.cc)
target_link_libraries(tests ${GTEST_LIBRARIES} gmock gmock_main pthread)
