
add_executable(main
	request.cpp
)

target_link_libraries(main
	cbd
)
target_include_directories(main
    PUBLIC
    "${PROJECT_SOURCE_DIR}/external/bm/src"
)