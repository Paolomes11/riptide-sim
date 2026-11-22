# io_root: the RDataFrame reader of the event files, the only part of the
# analysis that needs ROOT. ROOT is found through CMAKE_PREFIX_PATH or PATH
# (the directory holding root-config); no install path is hard-coded here.
find_package(ROOT 6.36 CONFIG REQUIRED COMPONENTS ROOTDataFrame)

add_library(riptide_io_root STATIC
    src/root_reader.cpp)
add_library(riptide::io_root ALIAS riptide_io_root)
target_include_directories(riptide_io_root PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_io_root
    PUBLIC riptide::io_core
    PRIVATE ROOT::ROOTDataFrame riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
