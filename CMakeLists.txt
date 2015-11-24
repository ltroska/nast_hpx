#################################
#Compile cdf_hpx
#@date 2015-11
#@author troska@ins.uni-bonn.de
#################################

#HPX requires this version
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

#C++-based project
project(cdf_hpx CXX)

find_package(HPX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/ /home/troska/opt/hpx_0.9.11/include 
			/home/troska/opt/boost_1_54_0/include
			/home/troska/opt/hwloc-1.11.1/include)

add_hpx_component(partition_server
    SOURCES src/grid/server/partition_server.cpp
    HEADERS src/grid/server/partition_server.hpp)

add_hpx_executable(debug
    ESSENTIAL
    SOURCES src/debug.cpp
    COMPONENT_DEPENDENCIES iostreams partition_server)
