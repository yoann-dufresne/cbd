

add_executable(unittests
	unittests.cpp
    Tests.cpp
)


target_link_libraries(unittests
	cbd
)


target_include_directories(unittests
    PUBLIC
    .
    "${PROJECT_SOURCE_DIR}/external/lest/include/lest"
    "${PROJECT_SOURCE_DIR}/external/bm/src"
)