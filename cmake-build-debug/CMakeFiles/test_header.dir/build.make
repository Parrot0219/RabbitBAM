# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/zhaozhan/CLionProjects/BamUniversalStatus

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/test_header.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/test_header.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_header.dir/flags.make

CMakeFiles/test_header.dir/test/test_header.cpp.o: CMakeFiles/test_header.dir/flags.make
CMakeFiles/test_header.dir/test/test_header.cpp.o: ../test/test_header.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_header.dir/test/test_header.cpp.o"
	/usr/local/Cellar/gcc@9/9.5.0/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_header.dir/test/test_header.cpp.o -c /Users/zhaozhan/CLionProjects/BamUniversalStatus/test/test_header.cpp

CMakeFiles/test_header.dir/test/test_header.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_header.dir/test/test_header.cpp.i"
	/usr/local/Cellar/gcc@9/9.5.0/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zhaozhan/CLionProjects/BamUniversalStatus/test/test_header.cpp > CMakeFiles/test_header.dir/test/test_header.cpp.i

CMakeFiles/test_header.dir/test/test_header.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_header.dir/test/test_header.cpp.s"
	/usr/local/Cellar/gcc@9/9.5.0/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zhaozhan/CLionProjects/BamUniversalStatus/test/test_header.cpp -o CMakeFiles/test_header.dir/test/test_header.cpp.s

# Object files for target test_header
test_header_OBJECTS = \
"CMakeFiles/test_header.dir/test/test_header.cpp.o"

# External object files for target test_header
test_header_EXTERNAL_OBJECTS =

test_header: CMakeFiles/test_header.dir/test/test_header.cpp.o
test_header: CMakeFiles/test_header.dir/build.make
test_header: CMakeFiles/test_header.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_header"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_header.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_header.dir/build: test_header
.PHONY : CMakeFiles/test_header.dir/build

CMakeFiles/test_header.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_header.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_header.dir/clean

CMakeFiles/test_header.dir/depend:
	cd /Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/zhaozhan/CLionProjects/BamUniversalStatus /Users/zhaozhan/CLionProjects/BamUniversalStatus /Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug /Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug /Users/zhaozhan/CLionProjects/BamUniversalStatus/cmake-build-debug/CMakeFiles/test_header.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_header.dir/depend

