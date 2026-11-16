add_executable(stats_tests
    numerics_test.cpp)
target_link_libraries(stats_tests PRIVATE
    riptide::stats riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(stats_tests PROPERTIES SKIP_RETURN_CODE 4)
