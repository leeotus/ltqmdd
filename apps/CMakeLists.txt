set(ORIGI_ALGO_NAME "orgnl")
set(MIXED_ALGO_NAME "mixed")
set(UPPER_ALGO_NAME "upper")
set(LOWER_ALGO_NAME "lower")
set(LTQMDDV1_TEST_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/inc")

add_executable("${ORIGI_ALGO_NAME}" orgnl-main.cpp)
add_executable("${MIXED_ALGO_NAME}" mixed-main.cpp)
add_executable("${UPPER_ALGO_NAME}" upper-main.cpp)
add_executable("${LOWER_ALGO_NAME}" lower-main.cpp)

target_link_libraries(
  ${ORIGI_ALGO_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings)
target_link_libraries(
  ${MIXED_ALGO_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings)
target_link_libraries(
  ${UPPER_ALGO_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings)
target_link_libraries(
  ${LOWER_ALGO_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings)

include_directories(src)
