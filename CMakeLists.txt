#################################
#Compile nast_hpx
#@date 2016-04
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(nast_hpx CXX)

find_package(HPX)

enable_testing()

#include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeSettings.txt)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/libs
			)

add_compile_options(-std=c++14)

add_library(pugixml ${CMAKE_CURRENT_SOURCE_DIR}/libs/pugixml/pugixml.cpp)
add_library(config src/io/config.cpp)
target_link_libraries(config pugixml)

# --------------- MAIN --------------- #

add_hpx_component(
    partition_server
    SOURCES src/grid/server/partition_server.cpp src/io/writer.cpp
    HEADERS src/grid/server/partition_server.hpp src/io/writer.hpp
)

add_hpx_component(
    stepper_server
    SOURCES src/stepper/server/stepper_server.cpp
    HEADERS src/stepper/server/stepper_server.hpp
    COMPILE_FLAGS -DCUSTOM_GRAIN_SIZE
    COMPONENT_DEPENDENCIES iostreams partition_server)

add_hpx_executable(
    main
    ESSENTIAL
    SOURCES src/main.cpp
    DEPENDENCIES config
    COMPONENT_DEPENDENCIES stepper_server
)

#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

#add_hpx_executable(
#    test_get_neighbor_cell
#    ESSENTIAL
#    SOURCES tests/test_get_neighbor_cell.cpp
#    COMPONENT_DEPENDENCIES partition_server
#)
#add_test(get_neighbor_cell test_get_neighbor_cell -t4)

