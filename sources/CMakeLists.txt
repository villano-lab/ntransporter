message("Building NT_Src")
add_executable(NT_Src NT_Src.cc)

target_include_directories(NT_Src PUBLIC "${CMAKE_SOURCE_DIR}") # for NTUtilities.hh

message("NT_Src built")