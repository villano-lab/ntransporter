message("Building NT_Src")
add_executable(NT_Src NT_Src.cc)

target_include_directories(NT_Src PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

message("NT_Src built")