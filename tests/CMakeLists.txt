

add_executable(unittests
	unittests.cpp
    #Tests.cpp
    serializeTest.cpp
)


target_link_libraries(unittests
	cbd
)


target_include_directories(unittests
    PUBLIC
    .
    "../external/lest/include/lest"
    "../external/bm/src"
)