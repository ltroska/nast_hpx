#################################
#Compile cdf_hpx
#@date 2015-11
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(nast_hpx CXX)

find_package(HPX)

enable_testing()

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeSettings.txt)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/libs)

add_library(pugixml ${CMAKE_CURRENT_SOURCE_DIR}/libs/pugixml/pugixml.cpp)

# --------------- MAIN --------------- #

add_hpx_component(
    partition_server
    SOURCES src/grid/server/partition_server.cpp
    HEADERS src/grid/server/partition_server.hpp
)

add_hpx_component(
    stepper_server
    SOURCES src/stepper/server/stepper_server.cpp
    HEADERS src/stepper/server/stepper_server.hpp
    COMPONENT_DEPENDENCIES iostreams partition_server)

add_hpx_executable(
    main
    ESSENTIAL
    SOURCES src/main.cpp src/io/manager.cpp src/solver/custom_chunk_solver.cpp
    DEPENDENCIES pugixml
    COMPONENT_DEPENDENCIES partition_server stepper_server
)

# --------------- UNIT TESTS --------------- #

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

# partition_data

add_hpx_executable(
    test_partition_data
    ESSENTIAL
    SOURCES tests/test_partition_data.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(partition_data_test test_partition_data)

# partition

add_hpx_executable(
    test_partition
    ESSENTIAL
    SOURCES tests/test_partition.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(partition_test test_partition)

