set(APP_TARGET_NAME "ltqmddv2")
set(LTQMDDV1_TEST_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/inc")

add_executable(${APP_TARGET_NAME} main.cpp)
target_link_libraries(
  ${APP_TARGET_NAME}  MQT::CoreDD MQT::CoreAlgorithms MQT::CoreCircuitOptimizer
                           MQT::ProjectOptions MQT::ProjectWarnings)

include_directories(src)