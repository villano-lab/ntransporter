message("Building PROCINFO")

add_executable(PROCINFO read_process_info.cc)
target_include_directories(PROCINFO PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh


target_link_libraries(PROCINFO PUBLIC "${Geant4_LIBRARIES}")
target_link_libraries(PROCINFO PUBLIC "${ROOT_LIBRARIES}")
target_link_libraries(PROCINFO SuperSim)

message("PROCINFO built")


