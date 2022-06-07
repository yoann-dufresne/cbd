cmake_minimum_required(VERSION 3.16)
project(testing)
add_compile_options(-Wall -Wextra)      #basic compiler options

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -Ofast -DNDEBUG") #Compilation flag

add_executable(main Tests.cpp)

add_subdirectory(cbl) #Where ConwayBromageLib contains CBl CMakeLists.txt 

target_link_libraries(main cbl)