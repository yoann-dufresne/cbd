cmake_minimum_required(VERSION 3.19)    #The minimal version of cmake we need


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



add_library(cbd ${SRCS} ${HEADERS})      #Use ConwayBromageLib as a library

target_include_directories(cbd PUBLIC .)

target_link_libraries(cbd
        sdsl
        divsufsort
        divsufsort64
)    #Specific flags from sdsl library

link_directories(cbd
        "${PROJECT_SOURCE_DIR}/external/sdsl-lite/lib"
        "${PROJECT_SOURCE_DIR}/external/sdsl-lite/external/libdivsufsort/lib"
)