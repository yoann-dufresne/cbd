cmake_minimum_required(VERSION 3.19)
project(ConwayBromageDatastructure)
add_compile_options(-Wall -Wextra)      #basic compiler options

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -Ofast -DNDEBUG") #Compilation flag


add_subdirectory(external)
add_subdirectory(cbd) #Where ConwayBromageLib contains CBl CMakeLists.txt 
add_subdirectory(apps)
add_subdirectory(tests)
