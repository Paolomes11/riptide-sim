add_executable(sim_validation_tests
    endf_format_test.cpp endf_reader_test.cpp check_stats_test.cpp
    expected_distributions_test.cpp sha256_test.cpp beam_geometry_test.cpp)
target_link_libraries(sim_validation_tests PRIVATE
    riptide::sim_validation riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_validation_tests PROPERTIES SKIP_RETURN_CODE 4)

# Reads the evaluated files under the data root; skipped (code 4) without them.
add_executable(sim_validation_reference_tests endf_reference_test.cpp)
target_link_libraries(sim_validation_reference_tests PRIVATE
    riptide::sim_validation riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_validation_reference_tests PROPERTIES SKIP_RETURN_CODE 4)
