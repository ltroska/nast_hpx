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

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeSettings.txt)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/libs
			)

add_compile_options(-std=c++14)

add_library(pugixml ${CMAKE_CURRENT_SOURCE_DIR}/libs/pugixml/pugixml.cpp)

# --------------- MAIN --------------- #

add_library(config src/io/config.cpp)
target_link_libraries(config pugixml)

add_hpx_component(
    partition_server
    SOURCES src/grid/server/partition_server.cpp
    HEADERS src/grid/server/partition_server.hpp
)

add_hpx_component(
	strategy
	SOURCES src/computation/custom_grain_size.cpp src/computation/with_for_each.cpp src/computation/cell_operations.cpp
	HEADERS src/computation/custom_grain_size.hpp src/computation/with_for_each.hpp src/computation/cell_operations.hpp
	COMPONENT_DEPENDENCIES partition_server
)

add_hpx_component(
    stepper_server
    SOURCES src/stepper/server/stepper_server.cpp
    HEADERS src/stepper/server/stepper_server.hpp
    COMPILE_FLAGS -DCUSTOM_GRAIN_SIZE
    COMPONENT_DEPENDENCIES iostreams partition_server strategy)

add_hpx_executable(
    main
    ESSENTIAL
    SOURCES src/main.cpp
    DEPENDENCIES config
    COMPONENT_DEPENDENCIES stepper_server
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

# partition_data

add_hpx_executable(
    test_partition_data
    ESSENTIAL
    SOURCES tests/test_partition_data.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(partition_data_test test_partition_data -t4)

# partition

add_hpx_executable(
    test_partition
    ESSENTIAL
    SOURCES tests/test_partition.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(partition_test test_partition -t4)

# get_neighbor_cell()

add_hpx_executable(
    test_get_neighbor_cell
    ESSENTIAL
    SOURCES tests/test_get_neighbor_cell.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(get_neighbor_cell test_get_neighbor_cell -t4)

# cgs::set_velocity_on_boundary()

add_hpx_executable(
    test_cgs_set_velocity_on_boundary
    ESSENTIAL
    SOURCES tests/test_cgs_set_velocity_on_boundary.cpp src/computation/custom_grain_size.cpp src/computation/cell_operations.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_set_velocity_on_boundary_test  test_cgs_set_velocity_on_boundary -t4)

