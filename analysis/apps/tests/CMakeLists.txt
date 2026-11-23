add_executable(analysis_cli_tests prepare_cli_test.cpp)
target_link_libraries(analysis_cli_tests PRIVATE
    riptide::analysis_cli riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(analysis_cli_tests PROPERTIES SKIP_RETURN_CODE 4)

add_executable(case_a_geant4_tests case_a_geant4_test.cpp)
target_link_libraries(case_a_geant4_tests PRIVATE
    riptide::case_a_geant4_core riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(case_a_geant4_tests PROPERTIES SKIP_RETURN_CODE 4)

# Smoke test of riptide_analysis_prepare on the runs written by sim_integration.
if(RIPTIDE_BUILD_ROOT_IO)
    add_test(NAME analysis_prepare_integration
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/run_prepare_it.sh"
            $<TARGET_FILE:riptide_analysis_prepare> "${PROJECT_SOURCE_DIR}")
    set_tests_properties(analysis_prepare_integration PROPERTIES
        LABELS integration
        FIXTURES_REQUIRED sim_integration
        SKIP_RETURN_CODE 4
        TIMEOUT 300)
    # Smoke test of riptide_case_a_geant4 on the same runs.
    add_test(NAME case_a_geant4_integration
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/run_case_a_it.sh"
            $<TARGET_FILE:riptide_case_a_geant4> $<TARGET_FILE:riptide_analysis_prepare>
            "${PROJECT_SOURCE_DIR}")
    set_tests_properties(case_a_geant4_integration PROPERTIES
        LABELS integration
        FIXTURES_REQUIRED sim_integration
        SKIP_RETURN_CODE 4
        TIMEOUT 300)
endif()
