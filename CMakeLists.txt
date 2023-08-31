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


set(CMAKE_CXX_STANDARD 14) # same as typically used for SuperSim

add_definitions(-DG4MULTITHREADED) # needed for correct run manager initialization

# for IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


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

# base supersim directory
set(SUPERSIM_DIR "/home/ajbiffl3/supersim" CACHE PATH "Path to base directory of SuperSim")

# SuperSim_include file
set(SUPERSIM_INCLUDE "SuperSim_include.txt" CACHE PATH "List of SuperSim directories to include")

message("Configuring SuperSim")
message(STATUS "SuperSim base directory: ${SUPERSIM_DIR}")
message(STATUS "Using SuperSim include list: ${SUPERSIM_INCLUDE}")


# read SuperSim_include into list
file(READ "${SUPERSIM_INCLUDE}" SUPERSIM_LIBRARIES_FILE)
string(REPLACE "\n" ";" SUPERSIM_LIBRARIES_FILE "${SUPERSIM_LIBRARIES_FILE}")
set(SUPERSIM_LIBRARIES "")
foreach (LIB ${SUPERSIM_LIBRARIES_FILE})
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
foreach(SUBDIR ${SUPERSIM_LIBRARIES})
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
		foreach(SOURCE_FILE ${SUBDIR_SOURCES})
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

# define SUPERSIM_LIB as a library with all source files and link to Geant and ROOT

set(SUPERSIM_INCLUDES "${SUPERSIM_DIR}")
foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	#target_include_directories(SUPERSIM_LIB PUBLIC "${SUPERSIM_DIR}/${SUBDIR}")
	LIST(APPEND SUPERSIM_INCLUDES "${SUPERSIM_DIR}/${SUBDIR}")
endforeach()


add_library(SUPERSIM::G4BASE SHARED IMPORTED)
add_library(SUPERSIM::APPS SHARED IMPORTED)
add_library(SUPERSIM::ACTIONS SHARED IMPORTED)
add_library(SUPERSIM::SOURCES SHARED IMPORTED)
add_library(SUPERSIM::INFO SHARED IMPORTED)
add_library(SUPERSIM::PHYSICS SHARED IMPORTED)
add_library(SUPERSIM::GEOMETRY SHARED IMPORTED)



#target_link_libraries(SUPERSIM_LIB PUBLIC ${Geant4_LIBRARIES})
#target_link_libraries(SUPERSIM_LIB PUBLIC ${ROOT_LIBRARIES})

#target_include_directories(SUPERSIM_LIB PUBLIC "${SUPERISM_DIR}")


set_target_properties(SUPERSIM::G4BASE PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSg4base.so"
					  )

set_target_properties(SUPERSIM::APPS PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSapps.so"
					  )

set_target_properties(SUPERSIM::ACTIONS PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSactions.so"
					  )

set_target_properties(SUPERSIM::SOURCES PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSsources.so"
					  )

set_target_properties(SUPERSIM::INFO PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSinfo.so"
					  )                                        

set_target_properties(SUPERSIM::PHYSICS PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSg4base.so"
					  )
                    
set_target_properties(SUPERSIM::GEOMETRY PROPERTIES 
					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSgeometry.so"
					  )


# set_target_properties(SUPERSIM_LIB PROPERTIES 
# 					  INTERFACE_INCLUDE_DIRECTORIES ${SUPERSIM_INCLUDES}
# 					  INTERFACE_SOURCES ${SUPERSOURCES}
# 					  IMPORTED_LOCATION "${SUPERSIM_DIR}/CDMSbuild/lib/Linux-g++/libCDMSapps.so"
# 					  )

add_library(SUPERSIM_LIB INTERFACE IMPORTED)
# set_target_properties(SUPERSIM_LIB PROPERTIES
#                       INTERFACE_LINK_LIBRARIES SUPERSIM::G4BASE SUPERSIM::APPS SUPERSIM::ACTIONS SUPERSIM::SOURCES SUPERSIM::INFO SUPERSIM::PHYSICS SUPERSIM::GEOMETRY
#                      )

set_property(TARGET SUPERSIM_LIB PROPERTY 
            INTERFACE_LINK_LIBRARIES
            SUPERSIM::G4BASE 
            SUPERSIM::APPS 
            SUPERSIM::ACTIONS 
            SUPERSIM::SOURCES 
            SUPERSIM::INFO 
            SUPERSIM::PHYSICS 
            SUPERSIM::GEOMETRY
            )


# ============================================
# add subdirectories
add_subdirectory(cross_sections)
add_subdirectory(sources)
add_subdirectory(boundary_conditions)





