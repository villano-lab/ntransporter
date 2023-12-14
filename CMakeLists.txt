# Base CMakeLists file - set up all executables

# Setup the project
cmake_minimum_required(VERSION 3.10)


set(CMAKE_BUILD_TYPE None)
message("cmake version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
message("cmake install prefix: ${CMAKE_INSTALL_PREFIX}")
message("source directory: ${CMAKE_CURRENT_SOURCE_DIR}")

project(NTransporter VERSION 1.0
        LANGUAGES CXX
        DESCRIPTION "Direct-coupled multigroup diffusion approximation with zeroth-order evaluated group constants")


set(CXX_STD 14 CACHE STRING "C++ standard to use")

set(CMAKE_CXX_STANDARD ${CXX_STD}) # 14 typically used for SuperSim, 17 used in conda installations of root

message("using C-std = C++${CMAKE_CXX_STANDARD}")

add_definitions(-DG4MULTITHREADED) # needed for correct run manager initialization

add_definitions(-DCDMS_USE_CVODE) # ensure CVODE packages are included

option(FORCE_LIBUUID_LINK "Force linking of libuuid shared object files with NT_XS and NT_XX" OFF)

option(IGNORE_WARNINGS "Ignore all warnings with -w flag" OFF)
if (IGNORE_WARNINGS)
	add_definitions(-w)
	message(WARNING "Ignoring all warnings")
endif()

# for IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


list(APPEND CMAKE_PREFIX_PATH "$ENV{G4INSTALL}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{CDMS_SUPERSIM}")

# =============================================
# Geant4

find_package(Geant4 REQUIRED)

if (Geant4_FOUND)
	message(STATUS "Libraries for Geant4 found")
        message(STATUS "Input directory: " ${Geant4_DIR})
	message(STATUS "Package Version: ${Geant4_VERSION}")
        message(STATUS "Include DIRS: ${Geant4_INCLUDE_DIRS}")
else()
	message(FATAL_ERROR "Failed to find Geant4.")
endif()

include("${Geant4_USE_FILE}")

# =============================================
# ROOT

find_package(ROOT REQUIRED)

if (ROOT_FOUND)
	message(STATUS "Libraries for ROOT found")
	message(STATUS "Input directory: ${ROOT_DIR}")
	message(STATUS "Package Version: ${ROOT_VERSION}")
	message(STATUS "Include DIRS: ${ROOT_INCLUDE_DIRS}")
else()
	message(FATAL_ERROR "Failed to find ROOT")
endif()

include("${ROOT_USE_FILE}")


# ============================================
# determine file extension for object/library files


if (${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")
	set(OBJ_EXT "dylib")
else()
	set(OBJ_EXT "so")
endif()
# ^ these conditions may not be perfect: can modify in future


# ============================================
# SuperSim


# SuperSim_include file
set(SUPERSIM_INCLUDE "SuperSim_include.txt" CACHE PATH "List of SuperSim directories to include")
set(SUPERSIM_EXCLUDE "SuperSim_exclude.txt" CACHE PATH "List of filenames to exclude from compilation")

message("Configuring SuperSim")
message(STATUS "SuperSim base directory: $ENV{CDMS_SUPERSIM}")
message(STATUS "Using SuperSim include list: ${SUPERSIM_INCLUDE}")
message(STATUS "Excluding ${SUPERSIM_EXCLUDE}")


# read SuperSim_include into list
file(READ "${SUPERSIM_INCLUDE}" SUPERSIM_LIBRARIES_FILE)
string(REPLACE "\n" ";" SUPERSIM_LIBRARIES_FILE "${SUPERSIM_LIBRARIES_FILE}")
set(SUPERSIM_LIBRARIES "")
foreach (LIB ${SUPERSIM_LIBRARIES_FILE})
	if ("${LIB}" MATCHES "^[^#]")
		list(APPEND SUPERSIM_LIBRARIES "${LIB}")
	endif()
endforeach()

message("looking for $ENV{G4WORKDIR}/lib/$ENV{G4SYSTEM}/*.${OBJ_EXT}")
#file(GLOB SUPERSIM_OBJ_ALL CONFIGURE_DEPENDS "$ENV{G4WORKDIR}/lib/$ENV{G4SYSTEM}/*.${OBJ_EXT}")


set(SUPERSIM_OBJ "")

#foreach(LIB "qhull" "qhull_p" "G4cmp" "CDMSsensitive" "CDMSfields" "CDMSreadout" "CDMSg4base" "CDMSsources" "CDMSgeometry" "CDMSfastdmc" "CDMSactions" "CDMShits" "CDMSinfo" "CDMSphysics" "CDMSapps")
foreach(LIB "CDMSapps" "CDMSphysics" "CDMSinfo" "CDMShits" "CDMSactions" "CDMSfastdmc" "CDMSgeometry" "CDMSsources" "CDMSg4base" "CDMSreadout" "CDMSfields" "CDMSsensitive" "G4cmp" "qhull_p" "qhull")
	file(GLOB FILE "$ENV{G4WORKDIR}/lib/$ENV{G4SYSTEM}/lib${LIB}.${OBJ_EXT}")
	list(APPEND SUPERSIM_OBJ ${FILE})
endforeach()


message("SUPERSIM_OBJ: ${SUPERSIM_OBJ}")

# execute script to configure SuperSim version header file in build directory
option(REGEN_CDMSVERSION "force re-generate CDMSVersion header in build directory" OFF)
if ((NOT EXISTS "${CMAKE_BINARY_DIR}/CDMSVersion.hh")
	OR (REGEN_CDMSVERSION))
	execute_process(
		COMMAND bash -c "source CDMSbuild/makeVersion-hh.sh CDMSbuild/CDMSVersion.template .cdms-gittag > ${CMAKE_BINARY_DIR}/CDMSVersion.hh"
		WORKING_DIRECTORY "$ENV{CDMS_SUPERSIM}"
		OUTPUT_VARIABLE SUPERSIM_VERSION_HH
		)
	message(STATUS "Configuring CDMSVersion.hh in ${CMAKE_BINARY_DIR}/CDMSVersion.hh")
endif()
set(REGEN_CDMSVERSION OFF CACHE INTERNAL "" FORCE)


# variable to hold list of filenames to skip
set(EXCLUDE "")

file(READ "${SUPERSIM_EXCLUDE}" SUPERSIM_EXCLUDE_FILE)
string(REPLACE "\n" ";" SUPERSIM_EXCLUDE_FILE "${SUPERSIM_EXCLUDE_FILE}")
foreach (FILE ${SUPERSIM_EXCLUDE_FILE})
	if ("${FILE}" MATCHES "^[^#]")
		list(APPEND EXCLUDE "${FILE}")
    endif()
endforeach()


# variable to hold list of SuperSim source files
set(SUPERSOURCES "")


option(EXCLUDE_VERBOSE "Print additional information about excluding files from linking/compilation" ON)

# read in all source files from directories listed in SuperSim_include, excluding files with names listed in SuperSim_exclude
foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	message(STATUS "Finding source in $ENV{CDMS_SUPERSIM}/${SUBDIR}")
	execute_process(
		COMMAND bash -c "find $ENV{CDMS_SUPERSIM}/${SUBDIR} -maxdepth 1 -name *.cc"
		WORKING_DIRECTORY $ENV{HOME}
		OUTPUT_VARIABLE SUBDIR_SOURCES
		)
	if ("${SUBDIR_SOURCES}" STREQUAL "")
		message("No source in $ENV{CDMS_SUPERSIM}/${SUBDIR}")
	else()
		string(REPLACE "\n" ";" SUBDIR_SOURCES "${SUBDIR_SOURCES}")
		foreach(SOURCE_FILE ${SUBDIR_SOURCES})
			#message("Source file: ${SOURCE_FILE}")
			set(SKIPTHIS FALSE)
			foreach(IGNORE ${EXCLUDE})
				#message("Checking to ignore ${IGNORE}")
				if ("${SOURCE_FILE}" MATCHES "/${IGNORE}$")
					#message("Matched ${IGNORE}")
					set(SKIPTHIS TRUE)
					break()
				endif()
			endforeach()
			#if (${SOURCE_FILE} IN_LIST EXCLUDE)
			#if (("${SUBDIR}" STREQUAL "CDMSapps") AND 
			#	(("${SOURCE_FILE}" MATCHES "CDMS_.*") OR 
			#	("${SOURCE_FILE}" MATCHES "readbackTiming.*") OR 
			#	("${SOURCE_FILE}" MATCHES "tetraNeighbors.*")))
			if (SKIPTHIS)
				if (EXCLUDE_VERBOSE)
					message("Skipping ${SOURCE_FILE}")
				endif()
			else()
				list(APPEND SUPERSOURCES "${SOURCE_FILE}")
			endif()
		endforeach()
	endif()
endforeach()

#file(READ "files.txt" SUPERSIM_FILES_FILE)
#string(REPLACE "\n" ";" SUPERSIM_FILES_FILE "${SUPERSIM_FILES_FILE}")
#foreach (FILE ${SUPERSIM_FILES_FILE})
#	if ("${FILE}" MATCHES "^[^#]")
#		list(APPEND SUPERSOURCES "${FILE}")
#	endif()
#endforeach()


# ============================================
# gather CVODE source files

#file(GLOB CVODE_LIBRARIES "$ENV{CVODE_HOME}/lib/*.dylib")

set(CVODE_LIB_FOLDER "lib" CACHE STRING "lib folder under CVODE_HOME containing CVODE objects (e.g., lib, lib64)")

message(STATUS "Linking CVODE objects in $ENV{CVODE_HOME}/${CVODE_LIB_FOLDER}")
#add_subdirectory("$ENV{CVODE_CMAKE_DIR}" build)
#add_library(CVODE OBJECT IMPORTED)

file(GLOB CVODE_OBJ CONFIGURE_DEPENDS "$ENV{CVODE_HOME}/${CVODE_LIB_FOLDER}/*.${OBJ_EXT}")


#set_target_properties(CVODE PROPERTIES IMPORTED_OBJECTS ${CVODE_OBJ})


#set(CVODE_FILES "cvode/cvode.c;cvode/cvode_direct.c;nvector/serial/nvector_serial.c;sunmatrix/dense/sunmatrix_dense.c;sunlinsol/dense/sunlinsol_dense.c")

#set(CVODE_SRC "")
#foreach(FILE ${CVODE_FILES})
#	message(STATUS "Adding CVODE files $ENV{CVODE_BASE}/${FILE}")
#	list(APPEND CVODE_SRC "$ENV{CVODE_BASE}/src/${FILE}")
#endforeach()

#add_library(CVODE SHARED "${CVODE_SRC}")
#add_library(CVODE SHARED "$ENV{CVODE_BASE}")
#add_subdirectory("$ENV{CVODE_BASE}/instdir/lib")
#add_library(CVODE OBJECT IMPORTED ${CVODE_SRC})
#set_target_properties(CVODE PROPERTIES IMPORTED_OBJECTS "$ENV{CVODE_BASE}/builddir/src")
#add_library(CVODE STATIC IMPORTED ${CVODE_SRC})
#set_target_properties(CVODE PROPERTIES IMPORTED_OBJECTS "$ENV{CVODE_BASE}/builddir/src")

# ============================================
# configure executables to all be put in base build directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# ============================================
# add subdirectories

option(SUBDIRS_VERBOSE "Print additional linking info for executables" OFF)


add_subdirectory(cross_sections)
add_subdirectory(sources)
add_subdirectory(boundary_conditions)


option(BUILD_PROCINFO "Build PROCINFO to print neutron physics process info" OFF)
if (BUILD_PROCINFO)
	add_subdirectory(process_reader)
endif()




