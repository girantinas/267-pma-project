# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /global/common/software/nersc/pm-2022q4/spack/linux-sles15-zen/cmake-3.24.3-k5msymx/bin/cmake

# The command to remove a file.
RM = /global/common/software/nersc/pm-2022q4/spack/linux-sles15-zen/cmake-3.24.3-k5msymx/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /global/homes/a/alehero/267-pma-project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /global/homes/a/alehero/267-pma-project/build

# Include any dependencies generated for this target.
include CMakeFiles/main_dist_pcsr.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main_dist_pcsr.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main_dist_pcsr.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main_dist_pcsr.dir/flags.make

CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o: CMakeFiles/main_dist_pcsr.dir/flags.make
CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o: /global/homes/a/alehero/267-pma-project/main_dist_pcsr.cpp
CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o: CMakeFiles/main_dist_pcsr.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/global/homes/a/alehero/267-pma-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o"
	/opt/cray/pe/craype/2.7.20/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o -MF CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o.d -o CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o -c /global/homes/a/alehero/267-pma-project/main_dist_pcsr.cpp

CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.i"
	/opt/cray/pe/craype/2.7.20/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /global/homes/a/alehero/267-pma-project/main_dist_pcsr.cpp > CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.i

CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.s"
	/opt/cray/pe/craype/2.7.20/bin/CC $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /global/homes/a/alehero/267-pma-project/main_dist_pcsr.cpp -o CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.s

# Object files for target main_dist_pcsr
main_dist_pcsr_OBJECTS = \
"CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o"

# External object files for target main_dist_pcsr
main_dist_pcsr_EXTERNAL_OBJECTS =

main_dist_pcsr: CMakeFiles/main_dist_pcsr.dir/main_dist_pcsr.cpp.o
main_dist_pcsr: CMakeFiles/main_dist_pcsr.dir/build.make
main_dist_pcsr: CMakeFiles/main_dist_pcsr.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/global/homes/a/alehero/267-pma-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main_dist_pcsr"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main_dist_pcsr.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main_dist_pcsr.dir/build: main_dist_pcsr
.PHONY : CMakeFiles/main_dist_pcsr.dir/build

CMakeFiles/main_dist_pcsr.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main_dist_pcsr.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main_dist_pcsr.dir/clean

CMakeFiles/main_dist_pcsr.dir/depend:
	cd /global/homes/a/alehero/267-pma-project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /global/homes/a/alehero/267-pma-project /global/homes/a/alehero/267-pma-project /global/homes/a/alehero/267-pma-project/build /global/homes/a/alehero/267-pma-project/build /global/homes/a/alehero/267-pma-project/build/CMakeFiles/main_dist_pcsr.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main_dist_pcsr.dir/depend
