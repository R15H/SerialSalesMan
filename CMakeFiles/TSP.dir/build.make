# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bruh/SalesmanSerial/SerialSalesMan

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bruh/SalesmanSerial/SerialSalesMan

# Include any dependencies generated for this target.
include CMakeFiles/TSP.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/TSP.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/TSP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TSP.dir/flags.make

CMakeFiles/TSP.dir/tsp.c.o: CMakeFiles/TSP.dir/flags.make
CMakeFiles/TSP.dir/tsp.c.o: tsp.c
CMakeFiles/TSP.dir/tsp.c.o: CMakeFiles/TSP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bruh/SalesmanSerial/SerialSalesMan/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/TSP.dir/tsp.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/TSP.dir/tsp.c.o -MF CMakeFiles/TSP.dir/tsp.c.o.d -o CMakeFiles/TSP.dir/tsp.c.o -c /home/bruh/SalesmanSerial/SerialSalesMan/tsp.c

CMakeFiles/TSP.dir/tsp.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TSP.dir/tsp.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/bruh/SalesmanSerial/SerialSalesMan/tsp.c > CMakeFiles/TSP.dir/tsp.c.i

CMakeFiles/TSP.dir/tsp.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TSP.dir/tsp.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/bruh/SalesmanSerial/SerialSalesMan/tsp.c -o CMakeFiles/TSP.dir/tsp.c.s

CMakeFiles/TSP.dir/queue.c.o: CMakeFiles/TSP.dir/flags.make
CMakeFiles/TSP.dir/queue.c.o: queue.c
CMakeFiles/TSP.dir/queue.c.o: CMakeFiles/TSP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bruh/SalesmanSerial/SerialSalesMan/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/TSP.dir/queue.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/TSP.dir/queue.c.o -MF CMakeFiles/TSP.dir/queue.c.o.d -o CMakeFiles/TSP.dir/queue.c.o -c /home/bruh/SalesmanSerial/SerialSalesMan/queue.c

CMakeFiles/TSP.dir/queue.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/TSP.dir/queue.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/bruh/SalesmanSerial/SerialSalesMan/queue.c > CMakeFiles/TSP.dir/queue.c.i

CMakeFiles/TSP.dir/queue.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/TSP.dir/queue.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/bruh/SalesmanSerial/SerialSalesMan/queue.c -o CMakeFiles/TSP.dir/queue.c.s

# Object files for target TSP
TSP_OBJECTS = \
"CMakeFiles/TSP.dir/tsp.c.o" \
"CMakeFiles/TSP.dir/queue.c.o"

# External object files for target TSP
TSP_EXTERNAL_OBJECTS =

bin/TSP: CMakeFiles/TSP.dir/tsp.c.o
bin/TSP: CMakeFiles/TSP.dir/queue.c.o
bin/TSP: CMakeFiles/TSP.dir/build.make
bin/TSP: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
bin/TSP: /usr/lib/x86_64-linux-gnu/libpthread.a
bin/TSP: CMakeFiles/TSP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bruh/SalesmanSerial/SerialSalesMan/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable bin/TSP"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TSP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TSP.dir/build: bin/TSP
.PHONY : CMakeFiles/TSP.dir/build

CMakeFiles/TSP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TSP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TSP.dir/clean

CMakeFiles/TSP.dir/depend:
	cd /home/bruh/SalesmanSerial/SerialSalesMan && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bruh/SalesmanSerial/SerialSalesMan /home/bruh/SalesmanSerial/SerialSalesMan /home/bruh/SalesmanSerial/SerialSalesMan /home/bruh/SalesmanSerial/SerialSalesMan /home/bruh/SalesmanSerial/SerialSalesMan/CMakeFiles/TSP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TSP.dir/depend

