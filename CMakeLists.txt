cmake_minimum_required(VERSION 3.0)

project(GAToolKit)

find_package(boost QUIET REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(GeneticsAlgorithem GeneticsAlgorithem.cpp)
add_executable(KnapsackGA KnapsackGA.cpp)
