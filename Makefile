# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_COMMAND = "/Applications/CMake 2.8-11.app/Contents/bin/cmake"

# The command to remove a file.
RM = "/Applications/CMake 2.8-11.app/Contents/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = "/Applications/CMake 2.8-11.app/Contents/bin/ccmake"

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bmetcalf/WorkSpace

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bmetcalf/WorkSpace

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	"/Applications/CMake 2.8-11.app/Contents/bin/ccmake" -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	"/Applications/CMake 2.8-11.app/Contents/bin/cmake" -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	cd /Users/bmetcalf/WorkSpace && $(CMAKE_COMMAND) -E cmake_progress_start /Users/bmetcalf/WorkSpace/CMakeFiles /Users/bmetcalf/WorkSpace/CosmoLib/CMakeFiles/progress.marks
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CMakeFiles/Makefile2 CosmoLib/all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/bmetcalf/WorkSpace/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CMakeFiles/Makefile2 CosmoLib/clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CMakeFiles/Makefile2 CosmoLib/preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CMakeFiles/Makefile2 CosmoLib/preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	cd /Users/bmetcalf/WorkSpace && $(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

# Convenience name for target.
CosmoLib/CMakeFiles/CosmoLib.dir/rule:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CMakeFiles/Makefile2 CosmoLib/CMakeFiles/CosmoLib.dir/rule
.PHONY : CosmoLib/CMakeFiles/CosmoLib.dir/rule

# Convenience name for target.
CosmoLib: CosmoLib/CMakeFiles/CosmoLib.dir/rule
.PHONY : CosmoLib

# fast build rule for target.
CosmoLib/fast:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/build
.PHONY : CosmoLib/fast

Cosmo/powerEH.o: Cosmo/powerEH.cpp.o
.PHONY : Cosmo/powerEH.o

# target to build an object file
Cosmo/powerEH.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEH.cpp.o
.PHONY : Cosmo/powerEH.cpp.o

Cosmo/powerEH.i: Cosmo/powerEH.cpp.i
.PHONY : Cosmo/powerEH.i

# target to preprocess a source file
Cosmo/powerEH.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEH.cpp.i
.PHONY : Cosmo/powerEH.cpp.i

Cosmo/powerEH.s: Cosmo/powerEH.cpp.s
.PHONY : Cosmo/powerEH.s

# target to generate assembly for a file
Cosmo/powerEH.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEH.cpp.s
.PHONY : Cosmo/powerEH.cpp.s

Cosmo/powerEHv2.o: Cosmo/powerEHv2.cpp.o
.PHONY : Cosmo/powerEHv2.o

# target to build an object file
Cosmo/powerEHv2.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEHv2.cpp.o
.PHONY : Cosmo/powerEHv2.cpp.o

Cosmo/powerEHv2.i: Cosmo/powerEHv2.cpp.i
.PHONY : Cosmo/powerEHv2.i

# target to preprocess a source file
Cosmo/powerEHv2.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEHv2.cpp.i
.PHONY : Cosmo/powerEHv2.cpp.i

Cosmo/powerEHv2.s: Cosmo/powerEHv2.cpp.s
.PHONY : Cosmo/powerEHv2.s

# target to generate assembly for a file
Cosmo/powerEHv2.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/Cosmo/powerEHv2.cpp.s
.PHONY : Cosmo/powerEHv2.cpp.s

cosmo.o: cosmo.cpp.o
.PHONY : cosmo.o

# target to build an object file
cosmo.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/cosmo.cpp.o
.PHONY : cosmo.cpp.o

cosmo.i: cosmo.cpp.i
.PHONY : cosmo.i

# target to preprocess a source file
cosmo.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/cosmo.cpp.i
.PHONY : cosmo.cpp.i

cosmo.s: cosmo.cpp.s
.PHONY : cosmo.s

# target to generate assembly for a file
cosmo.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/cosmo.cpp.s
.PHONY : cosmo.cpp.s

halo.o: halo.cpp.o
.PHONY : halo.o

# target to build an object file
halo.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/halo.cpp.o
.PHONY : halo.cpp.o

halo.i: halo.cpp.i
.PHONY : halo.i

# target to preprocess a source file
halo.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/halo.cpp.i
.PHONY : halo.cpp.i

halo.s: halo.cpp.s
.PHONY : halo.s

# target to generate assembly for a file
halo.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/halo.cpp.s
.PHONY : halo.cpp.s

nfw.o: nfw.cpp.o
.PHONY : nfw.o

# target to build an object file
nfw.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/nfw.cpp.o
.PHONY : nfw.cpp.o

nfw.i: nfw.cpp.i
.PHONY : nfw.i

# target to preprocess a source file
nfw.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/nfw.cpp.i
.PHONY : nfw.cpp.i

nfw.s: nfw.cpp.s
.PHONY : nfw.s

# target to generate assembly for a file
nfw.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/nfw.cpp.s
.PHONY : nfw.cpp.s

powerCDM.o: powerCDM.cpp.o
.PHONY : powerCDM.o

# target to build an object file
powerCDM.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/powerCDM.cpp.o
.PHONY : powerCDM.cpp.o

powerCDM.i: powerCDM.cpp.i
.PHONY : powerCDM.i

# target to preprocess a source file
powerCDM.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/powerCDM.cpp.i
.PHONY : powerCDM.cpp.i

powerCDM.s: powerCDM.cpp.s
.PHONY : powerCDM.s

# target to generate assembly for a file
powerCDM.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/powerCDM.cpp.s
.PHONY : powerCDM.cpp.s

utilities.o: utilities.cpp.o
.PHONY : utilities.o

# target to build an object file
utilities.cpp.o:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/utilities.cpp.o
.PHONY : utilities.cpp.o

utilities.i: utilities.cpp.i
.PHONY : utilities.i

# target to preprocess a source file
utilities.cpp.i:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/utilities.cpp.i
.PHONY : utilities.cpp.i

utilities.s: utilities.cpp.s
.PHONY : utilities.s

# target to generate assembly for a file
utilities.cpp.s:
	cd /Users/bmetcalf/WorkSpace && $(MAKE) -f CosmoLib/CMakeFiles/CosmoLib.dir/build.make CosmoLib/CMakeFiles/CosmoLib.dir/utilities.cpp.s
.PHONY : utilities.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... CosmoLib"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... Cosmo/powerEH.o"
	@echo "... Cosmo/powerEH.i"
	@echo "... Cosmo/powerEH.s"
	@echo "... Cosmo/powerEHv2.o"
	@echo "... Cosmo/powerEHv2.i"
	@echo "... Cosmo/powerEHv2.s"
	@echo "... cosmo.o"
	@echo "... cosmo.i"
	@echo "... cosmo.s"
	@echo "... halo.o"
	@echo "... halo.i"
	@echo "... halo.s"
	@echo "... nfw.o"
	@echo "... nfw.i"
	@echo "... nfw.s"
	@echo "... powerCDM.o"
	@echo "... powerCDM.i"
	@echo "... powerCDM.s"
	@echo "... utilities.o"
	@echo "... utilities.i"
	@echo "... utilities.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	cd /Users/bmetcalf/WorkSpace && $(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

