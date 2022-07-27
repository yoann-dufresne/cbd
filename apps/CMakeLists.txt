
add_executable(main
	main.cpp
)
add_executable(timerequest
    request.cpp
)
add_executable(nextrequest
    requestnext.cpp
)
add_executable(serializer
    serializer.cpp
)
add_executable(deserializertest
    deserializertest.cpp
)
target_link_libraries(main
	cbd
)
target_link_libraries(timerequest
	cbd
)
target_link_libraries(nextrequest
	cbd
)
target_link_libraries(serializer
	cbd
)
target_link_libraries(deserializertest
	cbd
)
set_target_properties( timerequest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties( nextrequest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties( serializer
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
set_target_properties(deserializertest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
target_include_directories(main
    PUBLIC
    "${PROJECT_SOURCE_DIR}/external/bm/src"
)