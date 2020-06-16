## Add sdsl-lite
set(sdsl_cmake_file ${CMAKE_CURRENT_LIST_DIR}/sdsl-lite/CMakeLists.txt)
if(NOT EXISTS ${sdsl_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/sdsl-lite
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
            OUTPUT_QUIET
            )
endif(NOT EXISTS ${sdsl_cmake_file})

add_subdirectory(sdsl-lite)
