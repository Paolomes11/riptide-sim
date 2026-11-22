add_executable(io_core_tests
    analysis_record_test.cpp
    response_test.cpp
    run_info_test.cpp
    selection_test.cpp)
target_link_libraries(io_core_tests PRIVATE
    riptide::io_core riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(io_core_tests PROPERTIES SKIP_RETURN_CODE 4)
