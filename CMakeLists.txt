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
	SOURCES src/computation/custom_grain_size.cpp
	HEADERS src/computation/custom_grain_size.hpp
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


