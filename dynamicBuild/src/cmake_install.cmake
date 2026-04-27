# Install script for directory: /home/planner/workshop_ljx/ltqmdd/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/ir/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/algorithms/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/circuit_optimizer/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/datastructures/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/dd/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/zx/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/na/cmake_install.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core" TYPE FILE FILES
    "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/mqt-core-config.cmake"
    "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/mqt-core-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/ir/libmqt-core-ir.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/algorithms/libmqt-core-algorithms.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/circuit_optimizer/libmqt-core-circuit-optimizer.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/datastructures/libmqt-core-ds.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/dd/libmqt-core-dd.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/zx/libmqt-core-zx.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/na/libmqt-core-na.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mqt-core" TYPE DIRECTORY FILES "/home/planner/workshop_ljx/ltqmdd/include/mqt-core/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmqt-core_Developmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core/mqt-core-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core/mqt-core-targets.cmake"
         "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/CMakeFiles/Export/share/cmake/mqt-core/mqt-core-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core/mqt-core-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core/mqt-core-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core" TYPE FILE FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/CMakeFiles/Export/share/cmake/mqt-core/mqt-core-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core" TYPE FILE FILES "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/src/CMakeFiles/Export/share/cmake/mqt-core/mqt-core-targets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/mqt-core" TYPE FILE FILES
    "/home/planner/workshop_ljx/ltqmdd/cmake/Cache.cmake"
    "/home/planner/workshop_ljx/ltqmdd/cmake/FindGMP.cmake"
    "/home/planner/workshop_ljx/ltqmdd/cmake/PackageAddTest.cmake"
    "/home/planner/workshop_ljx/ltqmdd/cmake/PreventInSourceBuilds.cmake"
    "/home/planner/workshop_ljx/ltqmdd/cmake/StandardProjectSettings.cmake"
    )
endif()

