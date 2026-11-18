add_executable(sim_tests
    units_consistency_test.cpp)
target_link_libraries(sim_tests PRIVATE
    riptide::common riptide::options riptide::warnings Catch2::Catch2WithMain ${Geant4_LIBRARIES})
catch_discover_tests(sim_tests PROPERTIES SKIP_RETURN_CODE 4)

# Integration test: real riptide_sim runs read back through the Geant4 ROOT
# reader. The reader is not discovered on its own: run_integration.sh makes the
# runs first and skips (code 4) when the data root (external SSD) is missing.
add_executable(sim_integration_reader
    integration_reader.cpp)
target_link_libraries(sim_integration_reader PRIVATE
    riptide::common riptide::sim_config riptide::sim_truth riptide::sim_reader riptide::options
    riptide::warnings Catch2::Catch2WithMain)
add_test(NAME sim_integration
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/run_integration.sh"
        $<TARGET_FILE:riptide_sim> $<TARGET_FILE:sim_integration_reader> "${PROJECT_SOURCE_DIR}")
set_tests_properties(sim_integration PROPERTIES
    LABELS integration
    SKIP_RETURN_CODE 4
    TIMEOUT 300)
