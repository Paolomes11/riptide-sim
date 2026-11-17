add_executable(sim_validation_tests
    endf_format_test.cpp endf_reader_test.cpp check_stats_test.cpp
    expected_distributions_test.cpp)
target_link_libraries(sim_validation_tests PRIVATE
    riptide::sim_validation riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_validation_tests PROPERTIES SKIP_RETURN_CODE 4)
