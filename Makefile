# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_COMMAND = /opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake

# The command to remove a file.
RM = /opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/runner/work/flip-fluid-simulation/flip-fluid-simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/runner/work/flip-fluid-simulation/flip-fluid-simulation/build

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target package
package: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Run CPack packaging tool..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cpack --config ./CPackConfig.cmake
.PHONY : package

# Special rule for the target package
package/fast: package
.PHONY : package/fast

# Special rule for the target package_source
package_source:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Run CPack packaging tool for source..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cpack --config ./CPackSourceConfig.cmake /home/runner/work/flip-fluid-simulation/flip-fluid-simulation/build/CPackSourceConfig.cmake
.PHONY : package_source

# Special rule for the target package_source
package_source/fast: package_source
.PHONY : package_source/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake cache editor..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Install the project..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Install the project..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing only the local directory..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing only the local directory..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing the project stripped..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing the project stripped..."
	/opt/hostedtoolcache/cmake/3.31.6/x64/cmake-3.31.6-linux-x86_64/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/runner/work/flip-fluid-simulation/flip-fluid-simulation/build/CMakeFiles /home/runner/work/flip-fluid-simulation/flip-fluid-simulation/build//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/runner/work/flip-fluid-simulation/flip-fluid-simulation/build/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named flip

# Build rule for target.
flip: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 flip
.PHONY : flip

# fast build rule for target.
flip/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/build
.PHONY : flip/fast

#=============================================================================
# Target rules for targets named uninstall

# Build rule for target.
uninstall: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 uninstall
.PHONY : uninstall

# fast build rule for target.
uninstall/fast:
	$(MAKE) $(MAKESILENT) -f vendor/raylib/CMakeFiles/uninstall.dir/build.make vendor/raylib/CMakeFiles/uninstall.dir/build
.PHONY : uninstall/fast

#=============================================================================
# Target rules for targets named raylib

# Build rule for target.
raylib: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 raylib
.PHONY : raylib

# fast build rule for target.
raylib/fast:
	$(MAKE) $(MAKESILENT) -f vendor/raylib/raylib/CMakeFiles/raylib.dir/build.make vendor/raylib/raylib/CMakeFiles/raylib.dir/build
.PHONY : raylib/fast

src/FlipFluid.o: src/FlipFluid.cpp.o
.PHONY : src/FlipFluid.o

# target to build an object file
src/FlipFluid.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/FlipFluid.cpp.o
.PHONY : src/FlipFluid.cpp.o

src/FlipFluid.i: src/FlipFluid.cpp.i
.PHONY : src/FlipFluid.i

# target to preprocess a source file
src/FlipFluid.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/FlipFluid.cpp.i
.PHONY : src/FlipFluid.cpp.i

src/FlipFluid.s: src/FlipFluid.cpp.s
.PHONY : src/FlipFluid.s

# target to generate assembly for a file
src/FlipFluid.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/FlipFluid.cpp.s
.PHONY : src/FlipFluid.cpp.s

src/main.o: src/main.cpp.o
.PHONY : src/main.o

# target to build an object file
src/main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/main.cpp.o
.PHONY : src/main.cpp.o

src/main.i: src/main.cpp.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/main.cpp.i
.PHONY : src/main.cpp.i

src/main.s: src/main.cpp.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/flip.dir/build.make CMakeFiles/flip.dir/src/main.cpp.s
.PHONY : src/main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... package"
	@echo "... package_source"
	@echo "... rebuild_cache"
	@echo "... uninstall"
	@echo "... flip"
	@echo "... raylib"
	@echo "... src/FlipFluid.o"
	@echo "... src/FlipFluid.i"
	@echo "... src/FlipFluid.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

