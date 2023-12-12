add_executable(NT_XS NT_XS.cc ${SUPERSOURCES})
add_executable(NT_XX NT_XX.cc ${SUPERSOURCES})

target_include_directories(NT_XS PUBLIC "$ENV{CVODE_HOME}/include")


#target_link_libraries(NT_XS PUBLIC CVODE)
#target_link_libraries(NT_XS PUBLIC "$ENV{CVODE_BASE}/instdir/lib")
#message(STATUS "NT_XS & NT_XX linking CVODE source in $ENV{CVODE_HOME}")
#message("CVODE_LIBRARIES = ${CVODE_LIBRARIES}")
#target_link_libraries(NT_XS PUBLIC ${CVODE_LIBRARIES})

target_include_directories(NT_XS PUBLIC "${CMAKE_SOURCE_DIR}") # for custom SuperSim_Main.hh
target_include_directories(NT_XX PUBLIC "${CMAKE_SOURCE_DIR}")

target_include_directories(NT_XS PUBLIC "${CMAKE_SOURCE_DIR}/include") # for NTUtilities.hh
target_include_directories(NT_XX PUBLIC "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(NT_XS PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC "${ROOT_LIBRARIES}")
target_link_libraries(NT_XX PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XX PUBLIC "${ROOT_LIBRARIES}")

# supersim
target_include_directories(NT_XS PUBLIC "$ENV{CDMS_SUPERSIM}")
target_include_directories(NT_XX PUBLIC "$ENV{CDMS_SUPERSIM}")

foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	if (SUBDIRS_VERBOSE)
		message(STATUS "NT_XS & NT_XX linking in $ENV{CDMS_SUPERSIM}/${SUBDIR}")
	endif()
	target_include_directories(NT_XS PUBLIC "$ENV{CDMS_SUPERSIM}/${SUBDIR}")
	target_include_directories(NT_XX PUBLIC "$ENV{CDMS_SUPERSIM}/${SUBDIR}")
endforeach()

target_include_directories(NT_XS PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh
target_include_directories(NT_XX PUBLIC "${CMAKE_BINARY_DIR}")

message("NT_XS & NT_XX linked")

# ============================================
# NT_DX
add_executable(NT_DX NT_DX.cc "${CMAKE_SOURCE_DIR}/src/MinimalUserInitialization.cc")

target_include_directories(NT_DX PUBLIC "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(NT_DX PUBLIC "${Geant4_LIBRARIES}")
#target_link_libraries(NT_DX PUBLIC "${CMAKE_SOURCE_DIR}/src")
