add_executable(PROCINFO read_process_info.cc ${SUPERSOURCES})
target_include_directories(PROCINFO PUBLIC "${CMAKE_SOURCE_DIR}") # for custom SuperSim_Main.hh
target_include_directories(PROCINFO PUBLIC "${CMAKE_SOURCE_DIR}/include") # for NTUtilities.hh


target_link_libraries(PROCINFO PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(PROCINFO PUBLIC "${ROOT_LIBRARIES}")

# supersim
target_include_directories(PROCINFO PUBLIC "$ENV{CDMS_SUPERSIM}")

target_include_directories(PROCINFO PUBLIC "$ENV{CVODE_HOME}/include")

target_link_libraries(PROCINFO PUBLIC ${CVODE_OBJ})

if (FORCE_LIBUUID_LINK)
	message(STATUS "Linking libuuid object as $ENV{LIBUUID_OBJ}")
	if (DEFINED ENV{LIBUUID_OBJ})
		target_link_libraries(PROCINFO PUBLIC "$ENV{LIBUUID_OBJ}")
	else()
		message(FATAL_ERROR "Error: env var LIBUUID_OBJ must be set with -DFORCE_LIBUUID_LINK=ON")
	endif()
endif()


foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	if (SUBDIRS_VERBOSE)
		message(STATUS "PROCINFO linking in $ENV{CDMS_SUPERSIM}/${SUBDIR}")
	endif()
	target_include_directories(PROCINFO PUBLIC "$ENV{CDMS_SUPERSIM}/${SUBDIR}")
endforeach()

target_include_directories(PROCINFO PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh

message("PROCINFO linked")


