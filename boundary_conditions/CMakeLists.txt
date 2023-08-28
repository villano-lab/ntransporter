message("Building NT_BC")
add_executable(NT_BC NT_BC.cc)

target_include_directories(NT_BC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}") # for NTUtilities.hh

message("NT_BC built")