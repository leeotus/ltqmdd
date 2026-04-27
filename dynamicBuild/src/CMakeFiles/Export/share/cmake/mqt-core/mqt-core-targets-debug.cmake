#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MQT::CoreIR" for configuration "Debug"
set_property(TARGET MQT::CoreIR APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreIR PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-ir.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreIR )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreIR "${_IMPORT_PREFIX}/lib/libmqt-core-ir.a" )

# Import target "MQT::CoreAlgorithms" for configuration "Debug"
set_property(TARGET MQT::CoreAlgorithms APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreAlgorithms PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-algorithms.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreAlgorithms )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreAlgorithms "${_IMPORT_PREFIX}/lib/libmqt-core-algorithms.a" )

# Import target "MQT::CoreCircuitOptimizer" for configuration "Debug"
set_property(TARGET MQT::CoreCircuitOptimizer APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreCircuitOptimizer PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-circuit-optimizer.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreCircuitOptimizer )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreCircuitOptimizer "${_IMPORT_PREFIX}/lib/libmqt-core-circuit-optimizer.a" )

# Import target "MQT::CoreDS" for configuration "Debug"
set_property(TARGET MQT::CoreDS APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreDS PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-ds.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreDS )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreDS "${_IMPORT_PREFIX}/lib/libmqt-core-ds.a" )

# Import target "MQT::CoreDD" for configuration "Debug"
set_property(TARGET MQT::CoreDD APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreDD PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-dd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreDD )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreDD "${_IMPORT_PREFIX}/lib/libmqt-core-dd.a" )

# Import target "MQT::CoreZX" for configuration "Debug"
set_property(TARGET MQT::CoreZX APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreZX PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-zx.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreZX )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreZX "${_IMPORT_PREFIX}/lib/libmqt-core-zx.a" )

# Import target "MQT::CoreNA" for configuration "Debug"
set_property(TARGET MQT::CoreNA APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(MQT::CoreNA PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libmqt-core-na.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MQT::CoreNA )
list(APPEND _IMPORT_CHECK_FILES_FOR_MQT::CoreNA "${_IMPORT_PREFIX}/lib/libmqt-core-na.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
