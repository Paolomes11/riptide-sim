# Integration test of the ROOT reader on the runs written by sim_integration
# (the fixture exists only with RIPTIDE_BUILD_SIM; without it the script skips
# when the runs are missing).
add_executable(io_root_integration
    root_reader_integration.cpp)
target_link_libraries(io_root_integration PRIVATE
    riptide::io_root riptide::io_core ROOT::ROOTDataFrame riptide::options riptide::warnings
    Catch2::Catch2WithMain)
add_test(NAME io_root_integration
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/run_io_root_it.sh"
        $<TARGET_FILE:io_root_integration> "${PROJECT_SOURCE_DIR}")
set_tests_properties(io_root_integration PROPERTIES
    LABELS integration
    FIXTURES_REQUIRED sim_integration
    SKIP_RETURN_CODE 4
    TIMEOUT 300)
