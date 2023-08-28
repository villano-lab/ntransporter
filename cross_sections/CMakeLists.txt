# Setup the project
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE None)
message("cmake install prefix: ${CMAKE_INSTALL_PREFIX}")
message("source directory: ${CMAKE_CURRENT_SOURCE_DIR}")

project(NT_XS) # the executable

add_definitions(-DG4MULTITHREADED) # needed for correct run manager initialization

list(APPEND CMAKE_PREFIX_PATH "$ENV{G4INSTALL}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
list(APPEND CMAKE_PREFIX_PATH "/home/ajbiffl3/supersim")

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
# SuperSim

message("Configuring SuperSim")

# base supersim directory
set(SUPERSIM_DIR "/home/ajbiffl3/supersim" CACHE PATH "Path to base directory of SuperSim")

# SuperSim_include file
set(SUPERSIM_INCLUDE "../SuperSim_include.txt" CACHE PATH "List of SuperSim directories to include")

# read SuperSim_include into list
file(READ "${SUPERSIM_INCLUDE}" SUPERSIM_LIBRARIES_FILE)
string(REPLACE "\n" ";" SUPERSIM_LIBRARIES_FILE "${SUPERSIM_LIBRARIES_FILE}")
set(SUPERSIM_LIBRARIES "")
foreach (LIB "${SUPERSIM_LIBRARIES_FILE}")
	if ("${LIB}" MATCHES "^[^#]")
		list(APPEND SUPERSIM_LIBRARIES "${LIB}")
	endif()
endforeach()


# execute script to configure SuperSim version header file in build directory
execute_process(
	COMMAND bash -c "source CDMSbuild/makeVersion-hh.sh CDMSbuild/CDMSVersion.template .cdms-gittag > ${CMAKE_BINARY_DIR}/CDMSVersion.hh"
	WORKING_DIRECTORY "${SUPERSIM_DIR}"
	OUTPUT_VARIABLE SUPERSIM_VERSION_HH
	)
message(STATUS "Configuring CDMSVersion.hh in ${CMAKE_BINARY_DIR}/CDMSVersion.hh")


# variable to hold list of SuperSim source files
set(SUPERSOURCES "")

# read in all source files from directories listed in SuperSim_include
foreach (SUBDIR "${SUPERSIM_LIBRARIES}")
	message(STATUS "Finding source in ${SUPERSIM_DIR}/${SUBDIR}")
	execute_process(
		COMMAND bash -c "find ${SUPERSIM_DIR}/${SUBDIR} -maxdepth 1 -name *.cc"
		WORKING_DIRECTORY $ENV{HOME}
		OUTPUT_VARIABLE SUBDIR_SOURCES
		)
	if ("${SUBDIR_SOURCES}" STREQUAL "")
		message("No source in ${SUPERSIM_DIR}/${SUBDIR}")
	else()
		string(REPLACE "\n" ";" SUBDIR_SOURCES "${SUBDIR_SOURCES}")
		foreach(SOURCE_FILE "${SUBDIR_SOURCES}")
			#message(${SOURCE_FILE})
			if (("${SUBDIR}" STREQUAL "CDMSapps") AND 
				(("${SOURCE_FILE}" MATCHES "CDMS_.*") OR 
				("${SOURCE_FILE}" MATCHES "readbackTiming.*") OR 
				("${SOURCE_FILE}" MATCHES "tetraNeighbors.*")))
				message("Skipping ${SOURCE_FILE}")
			else()
				list(APPEND SUPERSOURCES "${SOURCE_FILE}")
			endif()
		endforeach()
	endif()
endforeach()

# ============================================

# link everything to executable
add_executable(NT_XS NT_XS.cc "${SUPERSOURCES}")

target_link_libraries(NT_XS PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC "${ROOT_LIBRARIES}")

# include headers
target_include_directories(NT_XS PUBLIC "${SUPERSIM_DIR}")

foreach(SUBDIR "${SUPERSIM_LIBRARIES}")
	message(STATUS "Including and linking in ${SUPERSIM_DIR}/${SUBDIR}")
	target_include_directories(NT_XS PUBLIC "${SUPERSIM_DIR}/${SUBDIR}")
endforeach()

target_include_directories(NT_XS PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh


# remove cached variables to restore "default" behavior
unset(SUPERSIM_DIR CACHE)
unset(SUPERSIM_DIR)
unset(SUPERSIM_INCLUDE CACHE)
unset(SUPERSIM_INCLUDE)

