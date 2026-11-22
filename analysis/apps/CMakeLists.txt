# The analysis applications. Their command lines are a pure library, built and
# tested in every preset; the executables read the event files and so are
# built only with the ROOT reader.
add_library(riptide_analysis_cli STATIC src/prepare_cli.cpp)
add_library(riptide::analysis_cli ALIAS riptide_analysis_cli)
target_include_directories(riptide_analysis_cli PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_analysis_cli
    PUBLIC riptide::io_core
    PRIVATE riptide::ext_lyra riptide::options riptide::warnings)

if(RIPTIDE_BUILD_ROOT_IO)
    add_executable(riptide_analysis_prepare src/riptide_analysis_prepare.cpp)
    target_link_libraries(riptide_analysis_prepare PRIVATE
        riptide::analysis_cli riptide::io_core riptide::io_root riptide::common
        riptide::options riptide::warnings)
endif()

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
