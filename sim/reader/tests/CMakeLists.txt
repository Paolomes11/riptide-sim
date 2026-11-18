add_executable(sim_reader_tests
    run_reader_test.cpp)
target_link_libraries(sim_reader_tests PRIVATE
    riptide::sim_reader riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(sim_reader_tests PROPERTIES SKIP_RETURN_CODE 4)
