# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/alexa/Bureau/ConwayBromageLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alexa/Bureau/ConwayBromageLib

# Include any dependencies generated for this target.
include CMakeFiles/ConwayBromageLib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ConwayBromageLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ConwayBromageLib.dir/flags.make

CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o: CMakeFiles/ConwayBromageLib.dir/flags.make
CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o: Tests.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexa/Bureau/ConwayBromageLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o -c /home/alexa/Bureau/ConwayBromageLib/Tests.cpp

CMakeFiles/ConwayBromageLib.dir/Tests.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ConwayBromageLib.dir/Tests.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexa/Bureau/ConwayBromageLib/Tests.cpp > CMakeFiles/ConwayBromageLib.dir/Tests.cpp.i

CMakeFiles/ConwayBromageLib.dir/Tests.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ConwayBromageLib.dir/Tests.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexa/Bureau/ConwayBromageLib/Tests.cpp -o CMakeFiles/ConwayBromageLib.dir/Tests.cpp.s

CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o: CMakeFiles/ConwayBromageLib.dir/flags.make
CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o: ConwayBromageLib.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexa/Bureau/ConwayBromageLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o -c /home/alexa/Bureau/ConwayBromageLib/ConwayBromageLib.cpp

CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexa/Bureau/ConwayBromageLib/ConwayBromageLib.cpp > CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.i

CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexa/Bureau/ConwayBromageLib/ConwayBromageLib.cpp -o CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.s

# Object files for target ConwayBromageLib
ConwayBromageLib_OBJECTS = \
"CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o" \
"CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o"

# External object files for target ConwayBromageLib
ConwayBromageLib_EXTERNAL_OBJECTS =

ConwayBromageLib: CMakeFiles/ConwayBromageLib.dir/Tests.cpp.o
ConwayBromageLib: CMakeFiles/ConwayBromageLib.dir/ConwayBromageLib.cpp.o
ConwayBromageLib: CMakeFiles/ConwayBromageLib.dir/build.make
ConwayBromageLib: external/sdsl-lite/lib/libsdsl.a
ConwayBromageLib: external/sdsl-lite/external/libdivsufsort/lib/libdivsufsort.a
ConwayBromageLib: external/sdsl-lite/external/libdivsufsort/lib/libdivsufsort64.a
ConwayBromageLib: CMakeFiles/ConwayBromageLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alexa/Bureau/ConwayBromageLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ConwayBromageLib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ConwayBromageLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ConwayBromageLib.dir/build: ConwayBromageLib

.PHONY : CMakeFiles/ConwayBromageLib.dir/build

CMakeFiles/ConwayBromageLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ConwayBromageLib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ConwayBromageLib.dir/clean

CMakeFiles/ConwayBromageLib.dir/depend:
	cd /home/alexa/Bureau/ConwayBromageLib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alexa/Bureau/ConwayBromageLib /home/alexa/Bureau/ConwayBromageLib /home/alexa/Bureau/ConwayBromageLib /home/alexa/Bureau/ConwayBromageLib /home/alexa/Bureau/ConwayBromageLib/CMakeFiles/ConwayBromageLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ConwayBromageLib.dir/depend
