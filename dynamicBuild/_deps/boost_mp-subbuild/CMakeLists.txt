# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.22.1)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(boost_mp-populate NONE)



include(ExternalProject)
ExternalProject_Add(boost_mp-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "https://github.com/boostorg/multiprecision/archive/refs/tags/Boost_1_84_0.tar.gz"
                    SOURCE_DIR          "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/_deps/boost_mp-src"
                    BINARY_DIR          "/home/planner/workshop_ljx/ltqmdd/dynamicBuild/_deps/boost_mp-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


