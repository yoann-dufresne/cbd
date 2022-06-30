
add_executable(main
	main.cpp
)
add_executable(timerequest
    request.cpp
)
target_link_libraries(main
	cbd
)
target_link_libraries(timerequest
	cbd
)
set_target_properties( timerequest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
target_include_directories(main
    PUBLIC
    "${PROJECT_SOURCE_DIR}/external/bm/src"
)