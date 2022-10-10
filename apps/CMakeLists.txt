add_executable(buildIndex
    buildIndex.cpp
)
add_executable(loadIndex
    loadIndex.cpp
)
add_executable(queryIndex
    queryIndex.cpp
)

target_link_libraries(buildIndex
	cbd
)
target_link_libraries(loadIndex
	cbd
)
target_link_libraries(queryIndex
	cbd
)

set_target_properties(buildIndex
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(loadIndex
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(queryIndex
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)