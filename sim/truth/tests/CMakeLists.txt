add_executable(sim_truth_tests
    deposit_record_test.cpp
    event_row_test.cpp
    proton_light_test.cpp
    truth_record_test.cpp)
target_link_libraries(sim_truth_tests PRIVATE
    riptide::sim_truth riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_truth_tests PROPERTIES SKIP_RETURN_CODE 4)
