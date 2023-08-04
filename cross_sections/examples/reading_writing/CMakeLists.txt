# Setup the project
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE None)
message("cmake install prefix: ${CMAKE_INSTALL_PREFIX}")
message("source directory: ${CMAKE_CURRENT_SOURCE_DIR}")

project(NT_XSReader)

# =============================================

list(APPEND CMAKE_PREFIX_PATH "$ENV{G4INSTALL}")
list(APPEND CMAKE_PREFIX_PATH "/home/ajbiffl3/supersim")

find_package(Geant4 REQUIRED)

if (Geant4_FOUND)
	message(STATUS "Libraries for Geant4 found")
        message(STATUS "Input directory: " ${Geant4_DIR})
	message(STATUS "Package Version: " ${Geant4_VERSION})
        message(STATUS "Include DIRS: " ${Geant4_INCLUDE_DIRS})
else()
	message(FATAL_ERROR "Failed to find Geant4.")
endif()

include(${Geant4_USE_FILE})


link_libraries(${Geant4_LIBRARIES})

# =============================================

message("Configuring SuperSim")

# base supersim directory
set(SUPERSIM_DIR "/home/ajbiffl3/supersim" CACHE PATH "Path to base directory of SuperSim")

#add_library(supersim SHARED IMPORTED GLOBAL)
#set_target_properties(supersim PROPERTIES IMPORTED_LOCATION ${SUPERSIM_DIR}/GNUmakefile)

#if (TARGET supersim)
#	message(STATUS "Linking SuperSim at " ${SUPERSIM_DIR})
#else()
#	message(FATAL_ERROR "Failed to add SuperSim libraries. SUPERSIM_DIR was set as " ${SUPERSIM_DIR})
#endif()

set(SUPERSIM_INCLUDE "SuperSim_include.txt" CACHE PATH "List of SuperSim directories to include")

file(READ ${SUPERSIM_INCLUDE} SUPERSIM_LIBRARIES_FILE)


string(REPLACE "\n" ";" SUPERSIM_LIBRARIES_FILE ${SUPERSIM_LIBRARIES_FILE})
set(SUPERSIM_LIBRARIES "")


foreach (LIB ${SUPERSIM_LIBRARIES_FILE})
	if (${LIB} MATCHES "^[^#]")
		list(APPEND SUPERSIM_LIBRARIES ${LIB})
	endif()
endforeach()


set(SUPERSOURCES "")


foreach (SUBDIR ${SUPERSIM_LIBRARIES})
	message(STATUS "Finding headers in " ${SUPERSIM_DIR}/${SUBDIR})
	#file(GLOB_RECURSE SUBDIR_HEADERS CONFIGURE_DEPENDS "${SUPERSIM_DIR}/${SUBDIR}/*.h?")
	execute_process(
		COMMAND bash -c "find ${SUPERSIM_DIR}/${SUBDIR} -name *.hh"
		WORKING_DIRECTORY $ENV{HOME}
		OUTPUT_VARIABLE SUBDIR_HEADERS
		) 
	string(REPLACE "\n" ";" SUBDIR_HEADERS ${SUBDIR_HEADERS})
	#message(STATUS "Found" ${SUBDIR_HEADERS})
	foreach(HEADER_FILE ${SUBDIR_HEADERS})
		#message("hfile:" ${HEADER_FILE})
		list(APPEND SUPERSOURCES ${HEADER_FILE})
	endforeach()
endforeach()

#message("SUPERSOURCES:\n" ${SUPERSOURCES})
#foreach(FILE ${SUPERSOURCES})
#	message(${FILE})
#endforeach()

#link_libraries(${Geant4_LIBRARIES})#, supersim)

add_executable(NT_XSReader read_xs.cc ${SUPERSOURCES})

foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	message(STATUS "Including and linking in " ${SUPERSIM_DIR}/${SUBDIR})
	target_include_directories(NT_XSReader PUBLIC ${SUPERSIM_DIR}/${SUBDIR})
endforeach()



# remove cached variables to restore "default" behavior
unset(SUPERSIM_DIR CACHE)
unset(SUPERSIM_DIR)
unset(SUPERSIM_INCLUDE CACHE)
unset(SUPERSIM_INCLUDE)

