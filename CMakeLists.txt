cmake_minimum_required(VERSION 3.0)

project(GAToolKit)

if (WIN32)
    find_package(boost QUIET REQUIRED)
endif(WIN32)
if (UNIX)
    find_package(Boost QUIET REQUIRED)
endif(UNIX)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(GeneticsAlgorithem GeneticsAlgorithem.cpp)
add_executable(KnapsackGA KnapsackGA.cpp)
add_executable(2DKnapsackGA 2DKnapsackGA.cpp)
