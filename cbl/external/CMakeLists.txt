## Add sdsl-lite
set(sdsl_cmake_file ${CMAKE_CURRENT_LIST_DIR}/sdsl-lite/CMakeLists.txt) #Variable
if(NOT EXISTS ${sdsl_cmake_file})   #Verify if it is possible de find the sdsl-lite CMake
    EXECUTE_PROCESS(COMMAND git submodule update --init external/sdsl-lite  #Verify updates of the submodule
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
            OUTPUT_QUIET
            )
endif(NOT EXISTS ${sdsl_cmake_file})

add_subdirectory(sdsl-lite) #Considered sdsl-lite as a subdirectory of our library

#Don't need it for lest because lest is just one .h
