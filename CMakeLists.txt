cmake_minimum_required(VERSION 3.12)    #The minimal version of cmake we need
#Here we need 3.16 or higher version
project(CBL)               #The name of the projet, it will be the executable name

add_compile_options(-Wall -Wextra)      #basic compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG") #Compilation flag
#Variables which represents compilation elements
#One for .cpp elements
set(SRCS
        Tests.cpp
        ConwayBromageLib.cpp)
#One for .h
set(HEADERS
        ConwayBromageLib.h)
link_directories(ConwayBromageLib ./sdsl-lite/lib ./sdsl-lite/external/libdivsufsort/lib)    #Where we put .cpp from sdsl library
add_executable(ConwayBromageLib ${SRCS} ${HEADERS}) #Create executable ConwayBromageLib
add_subdirectory(external)      #Add subdirectories like sdsl and lest, call the CMakeLists.txt from external directory
target_include_directories(ConwayBromageLib PUBLIC ./external/lest/include)   #Where we put .h from lest libraries : It needs to be after add_executable !!!
target_link_libraries(ConwayBromageLib sdsl divsufsort divsufsort64)    #Specific flags from sdsl library

#Original compilation line without using CMake :
#g++ -Wall -Wextra -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG -I ~/include -L ~/lib -o Tests.exe Tests.cpp ConwayBromageLib.cpp -lsdsl -ldivsufsort -ldivsufsort64
