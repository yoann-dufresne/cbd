cmake_minimum_required(VERSION 3.16)
project(ConwayBromageLib)

set(CMAKE_CXX_STANDARD 11)
include_directories("~/include/")
add_executable(ConwayBromageLib main.cpp ConwayBromageLib.cpp ConwayBromageLib.h)
target_link_directories(ConwayBromageLib PRIVATE ~/lib/)