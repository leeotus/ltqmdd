set(NODEPAR_UNITEST_NAME "mnodepar-test")

add_executable("${NODEPAR_UNITEST_NAME}" test_mnodepar.cpp)

target_link_libraries(
  ${NODEPAR_UNITEST_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings gtest_main gmock_main)
