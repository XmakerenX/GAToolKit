cmake_minimum_required(VERSION 3.0)

project(gatoolkit)

add_executable(gatoolkit main.cpp)

install(TARGETS gatoolkit RUNTIME DESTINATION bin)
