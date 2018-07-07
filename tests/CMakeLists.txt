FUNCTION(add_test_executable TEST_NAME)
    add_executable(${TEST_NAME} ${ARGN})
    target_compile_features(${TEST_NAME} PRIVATE cxx_range_for)
    add_test(${TEST_NAME} ${TEST_NAME})
ENDFUNCTION()

include_directories(${CMAKE_CURRENT_DIRECTORY})
include_directories("${PROJECT_SOURCE_DIR}")

enable_testing()

add_test_executable(test_wavefront     wavefront.cpp)
add_test_executable(test_wavefront_bvh wavefront_bvh.cpp)
add_test_executable(test_plate         plate.cpp)
add_test_executable(test_bounding_box  bounding_box.cpp)
