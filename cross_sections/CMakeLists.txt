message("Building NT_XS")

add_executable(NT_XS NT_XS.cc ${SUPERSOURCES})
target_include_directories(NT_XS PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh


target_link_libraries(NT_XS PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC "${ROOT_LIBRARIES}")

# supersim
target_include_directories(NT_XS PUBLIC "$ENV{CDMS_SUPERSIM}")

foreach(SUBDIR ${SUPERSIM_LIBRARIES})
	message(STATUS "NT_XS linking in $ENV{CDMS_SUPERSIM}/${SUBDIR}")
	target_include_directories(NT_XS PUBLIC "$ENV{CDMS_SUPERSIM}/${SUBDIR}")
endforeach()

target_include_directories(NT_XS PUBLIC "${CMAKE_BINARY_DIR}") # for CDMSVersion.hh

message("NT_XS built")


