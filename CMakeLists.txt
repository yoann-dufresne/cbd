cmake_minimum_required(VERSION 3.12)    #The minimal version of cmake we need
#Here we need 3.12 or higher version
project(ConwayBromageLib
        VERSION
        DESCRIPTION "ConwayBromageLib : an implementation of Conway and Bromage succinct data structure for assembling large genomes"
        HOMEPAGE_URL "https://github.com/yoann-dufresne/ConwayBromageLib.git"
        LANGUAGES CXX)          #information about the project
#Variables which represents compilation elements
#One for .cpp elements
set(SRCS
        ConwayBromageLib.cpp)
#One for .h
set(HEADERS
        ConwayBromageLib.h)
link_directories(cbl ./sdsl-lite/lib ./sdsl-lite/external/libdivsufsort/lib)    #Where we put .cpp from sdsl library
add_library(cbl ${SRCS} ${HEADERS})      #Use ConwayBromageLib as a library
target_include_directories(cbl PUBLIC .)
add_subdirectory(external)
target_link_libraries(cbl sdsl divsufsort divsufsort64)    #Specific flags from sdsl library
