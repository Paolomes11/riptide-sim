# io_core: everything between the event files and the analysis that does not
# need ROOT (run.json, selection, HP conversion, resolution, analysis.json),
# so that it builds and is tested in CI too.
add_library(riptide_io_core STATIC
    src/analysis_record.cpp
    src/event_columns.cpp
    src/response.cpp
    src/response_summary.cpp
    src/run_info.cpp
    src/selection.cpp)
add_library(riptide::io_core ALIAS riptide_io_core)
target_include_directories(riptide_io_core PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_io_core
    PUBLIC riptide::common riptide::ext_json
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
