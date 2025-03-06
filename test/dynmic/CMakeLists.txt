set(NODEPAR_UNITEST_NAME "mnodepar-test")
set(SIFTING_ALGO_UNITEST_NAME "siftingalgo-test")

add_executable("${NODEPAR_UNITEST_NAME}" test_mnodepar.cpp)
add_executable("${SIFTING_ALGO_UNITEST_NAME}" test_sifting_algo.cpp)

target_link_libraries(
  ${NODEPAR_UNITEST_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings gtest_main gmock_main)

target_link_libraries(
  ${SIFTING_ALGO_UNITEST_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings gtest_main gmock_main)
