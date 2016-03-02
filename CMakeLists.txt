#################################
#Compile nast_hpx
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
    SOURCES src/main.cpp src/io/manager.cpp src/computation/custom_grain_size.cpp
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
    SOURCES tests/test_cgs_set_velocity_on_boundary.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_set_velocity_on_boundary_test  test_cgs_set_velocity_on_boundary -t4)

# cgs::set_pressure_on_boundary()

add_hpx_executable(
    test_cgs_set_pressure_on_boundary
    ESSENTIAL
    SOURCES tests/test_cgs_set_pressure_on_boundary.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_set_pressure_on_boundary_test  test_cgs_set_pressure_on_boundary -t4)

# cgs::compute_fg()

add_hpx_executable(
    test_cgs_compute_fg
    ESSENTIAL
    SOURCES tests/test_cgs_compute_fg.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_compute_fg_test  test_cgs_compute_fg -t4)

# cgs::compute_rhs()

add_hpx_executable(
    test_cgs_compute_rhs
    ESSENTIAL
    SOURCES tests/test_cgs_compute_rhs.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_compute_rhs_test  test_cgs_compute_rhs -t4)

# cgs::sor_cycle()

add_hpx_executable(
    test_cgs_sor_cycle
    ESSENTIAL
    SOURCES tests/test_cgs_sor_cycle.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_sor_cycle_test  test_cgs_sor_cycle -t4)

# cgs::update_velocities()

add_hpx_executable(
    test_cgs_update_velocities
    ESSENTIAL
    SOURCES tests/test_cgs_update_velocities.cpp src/computation/custom_grain_size.cpp
    COMPONENT_DEPENDENCIES partition_server
)
add_test(cgs_update_velocities_test  test_cgs_update_velocities -t4)

