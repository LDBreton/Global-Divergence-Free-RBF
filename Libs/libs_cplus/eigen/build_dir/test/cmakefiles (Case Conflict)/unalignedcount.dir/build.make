# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir"

# Include any dependencies generated for this target.
include test/CMakeFiles/unalignedcount.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/unalignedcount.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/unalignedcount.dir/flags.make

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o: test/CMakeFiles/unalignedcount.dir/flags.make
test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o: ../test/unalignedcount.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o"
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o -c "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/test/unalignedcount.cpp"

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/unalignedcount.dir/unalignedcount.cpp.i"
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/test/unalignedcount.cpp" > CMakeFiles/unalignedcount.dir/unalignedcount.cpp.i

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/unalignedcount.dir/unalignedcount.cpp.s"
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/test/unalignedcount.cpp" -o CMakeFiles/unalignedcount.dir/unalignedcount.cpp.s

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.requires:

.PHONY : test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.requires

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.provides: test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/unalignedcount.dir/build.make test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.provides.build
.PHONY : test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.provides

test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.provides.build: test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o


# Object files for target unalignedcount
unalignedcount_OBJECTS = \
"CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o"

# External object files for target unalignedcount
unalignedcount_EXTERNAL_OBJECTS =

test/unalignedcount: test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o
test/unalignedcount: test/CMakeFiles/unalignedcount.dir/build.make
test/unalignedcount: test/CMakeFiles/unalignedcount.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable unalignedcount"
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/unalignedcount.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/unalignedcount.dir/build: test/unalignedcount

.PHONY : test/CMakeFiles/unalignedcount.dir/build

test/CMakeFiles/unalignedcount.dir/requires: test/CMakeFiles/unalignedcount.dir/unalignedcount.cpp.o.requires

.PHONY : test/CMakeFiles/unalignedcount.dir/requires

test/CMakeFiles/unalignedcount.dir/clean:
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" && $(CMAKE_COMMAND) -P CMakeFiles/unalignedcount.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/unalignedcount.dir/clean

test/CMakeFiles/unalignedcount.dir/depend:
	cd "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb" "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/test" "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir" "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test" "/home/porufes/Dropbox/Programas Numericos/FBR C++/eigen-eigen-5a0156e40feb/build_dir/test/CMakeFiles/unalignedcount.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : test/CMakeFiles/unalignedcount.dir/depend

