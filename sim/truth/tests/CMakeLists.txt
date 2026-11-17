add_executable(sim_truth_tests
    truth_record_test.cpp)
target_link_libraries(sim_truth_tests PRIVATE
    riptide::sim_truth riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_truth_tests PROPERTIES SKIP_RETURN_CODE 4)
