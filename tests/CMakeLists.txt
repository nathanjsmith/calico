FUNCTION(add_test_executable TEST_NAME)
    add_executable(${TEST_NAME} ${ARGN})
    target_compile_features(${TEST_NAME} PRIVATE cxx_range_for)
    target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(${TEST_NAME} ${TEST_NAME})
ENDFUNCTION()

find_package(Boost COMPONENTS system unit_test_framework REQUIRED)

include_directories(${CMAKE_CURRENT_DIRECTORY})
include_directories("${PROJECT_SOURCE_DIR}")

enable_testing()

add_test_executable(test_wavefront     wavefront.cpp)
add_test_executable(test_plate         plate.cpp)
add_test_executable(test_bounding_box  bounding_box.cpp)
