cmake_minimum_required(VERSION 3.19)
project(ConwayBromageDatastructure)
add_compile_options(-Wall -Wextra)      #basic compiler options

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++17 -DNDEBUG -Ofast ") #Compilation flag


add_subdirectory(external)
add_subdirectory(cbd) #Where ConwayBromageLib contains CBl CMakeLists.txt 
add_subdirectory(apps)
add_subdirectory(tests)
