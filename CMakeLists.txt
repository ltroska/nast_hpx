#################################
#Compile cdf_hpx
#@date 2015-11
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(cfd_hpx CXX)

find_package(HPX)

include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeSettings.txt)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/libs)

add_library(pugixml ${CMAKE_CURRENT_SOURCE_DIR}/libs/pugixml/pugixml.cpp)

add_hpx_component(partition_server
    SOURCES src/grid/server/partition_server.cpp
    HEADERS src/grid/server/partition_server.hpp)

add_hpx_executable(debug
    ESSENTIAL
    SOURCES src/debug.cpp src/io/config_reader.cpp
    DEPENDENCIES pugixml
    COMPONENT_DEPENDENCIES iostreams partition_server)

