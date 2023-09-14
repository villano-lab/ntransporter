message("Building NT_XS")

add_executable(NT_XS NT_XS.cc)
target_include_directories(NT_XS PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh


target_link_libraries(NT_XS PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC "${ROOT_LIBRARIES}")
target_link_libraries(NT_XS PUBLIC SuperSim)

message("NT_XS built")


