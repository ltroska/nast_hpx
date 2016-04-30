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

#enable_testing()

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeSettings.txt)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/libs)

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
	cgs_strategy
	SOURCES src/computation/custom_grain_size.cpp src/computation/cell_operations.cpp
	HEADERS src/computation/custom_grain_size.hpp src/computation/cell_operations.hpp
	COMPONENT_DEPENDENCIES partition_server
)

add_hpx_component(
	wfe_strategy
	SOURCES src/computation/with_for_each.cpp src/computation/cell_operations.cpp
	HEADERS src/computation/with_for_each.hpp src/computation/cell_operations.hpp
	COMPONENT_DEPENDENCIES partition_server
)

add_hpx_component(
    stepper_server_cgs
    SOURCES src/stepper/server/stepper_server.cpp
    HEADERS src/stepper/server/stepper_server.hpp
    COMPILE_FLAGS -DCUSTOM_GRAIN_SIZE
    COMPONENT_DEPENDENCIES iostreams partition_server cgs_strategy)

add_hpx_component(
    stepper_server_wfe
    SOURCES src/stepper/server/stepper_server.cpp
    HEADERS src/stepper/server/stepper_server.hpp
    COMPILE_FLAGS
    COMPONENT_DEPENDENCIES iostreams partition_server wfe_strategy)

add_hpx_executable(
    cgs
    ESSENTIAL
    SOURCES src/main.cpp
    DEPENDENCIES config
    COMPONENT_DEPENDENCIES stepper_server_cgs
)

add_hpx_executable(
    wfe
    ESSENTIAL
    SOURCES src/main.cpp
    DEPENDENCIES config
    COMPONENT_DEPENDENCIES stepper_server_wfe
)
