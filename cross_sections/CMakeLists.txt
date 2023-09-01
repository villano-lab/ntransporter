message("Building NT_XS")

#add_executable(NT_XS NT_XS.cc ${SUPERSOURCES})
add_executable(NT_XS NT_XS.cc)
target_include_directories(NT_XS PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh


target_link_libraries(NT_XS PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC "${ROOT_LIBRARIES}")

# supersim
#message("Linking ${SUPERSIM_LIB_LIBRARIES}")
add_dependencies(NT_XS SUPERSIM_G4BASE SUPERSIM_APPS SUPERSIM_ACTIONS SUPERSIM_SOURCES SUPERSIM_INFO SUPERSIM_PHYSICS SUPERSIM_GEOMETRY)
target_link_libraries(NT_XS PUBLIC SUPERSIM_LIB)
target_include_directories(NT_XS PUBLIC "${SUPERSIM_DIR}")

foreach(SUBDIR ${SUPERSIM_LIBRARIES})
 	message(STATUS "NT_XS linking in ${SUPERSIM_DIR}/${SUBDIR}")
 	target_include_directories(NT_XS PUBLIC "${SUPERSIM_DIR}/${SUBDIR}")
endforeach()

target_include_directories(NT_XS PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh

message("NT_XS built")


