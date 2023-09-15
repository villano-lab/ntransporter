message("Building PROCINFO")

add_executable(PROCINFO read_process_info.cc ${SUPERSOURCES})
target_include_directories(PROCINFO PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh


target_link_libraries(PROCINFO PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(PROCINFO PUBLIC "${ROOT_LIBRARIES}")

# supersim
target_include_directories(PROCINFO PUBLIC "${SUPERSIM_DIR}")
target_include_directories(PROCINFO PUBLIC "${UUID_INCLUDEDIR}")

foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	message(STATUS "PROCINFO linking in ${SUPERSIM_DIR}/${SUBDIR}")
	target_include_directories(PROCINFO PUBLIC "${SUPERSIM_DIR}/${SUBDIR}")
endforeach()

target_include_directories(PROCINFO PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh

message("PROCINFO built")


