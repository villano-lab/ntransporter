add_executable(NT_BC NT_BC.cc)

target_include_directories(NT_BC PUBLIC "${CMAKE_SOURCE_DIR}/include") # for NTUtilities.hh

message("NT_BC linked")