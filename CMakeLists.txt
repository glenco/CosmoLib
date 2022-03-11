cmake_minimum_required(VERSION 3.21)

####
# project settings
####

project(CosmoLib)

set(CosmoLib_VERSION_MAJOR 1)
set(CosmoLib_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
	include/cosmo.h
	include/halo.h
	include/powerCDMHM.h
	include/utilities.h
)


####
# sources
####

set(SOURCES
	Cosmo/powerEH.cpp
	Cosmo/powerEHv2.cpp
	cosmo.cpp
	halo.cpp
	nfw.cpp
	powerCDM.cpp
	powerCDMHM.cpp
	utilities.cpp
)


####
# dependencies
####

if(ENABLE_GSL)
	include_directories(${GSL_INCLUDE_DIRS})
endif()

include_directories(${NR_SOURCE_DIR}/include)
include_directories(include)


####
# build
####

add_library(CosmoLib ${HEADERS} ${SOURCES})

add_dependencies(CosmoLib NR)

