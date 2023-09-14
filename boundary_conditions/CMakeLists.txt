message("Building NT_BC")
add_executable(NT_BC NT_BC.cc)

target_include_directories(NT_BC PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

message("NT_BC built")