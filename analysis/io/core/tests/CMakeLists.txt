add_executable(io_core_tests
    analysis_record_test.cpp
    event_columns_test.cpp
    response_summary_test.cpp
    response_test.cpp
    run_info_test.cpp
    selection_test.cpp)
target_link_libraries(io_core_tests PRIVATE
    riptide::io_core riptide::options riptide::warnings Catch2::Catch2WithMain)
# The shipped configs are read as they are, to test the run.json source reader
# against the files a real run copies.
target_compile_definitions(io_core_tests PRIVATE
    RIPTIDE_CONFIG_DIR="${PROJECT_SOURCE_DIR}/config")
catch_discover_tests(io_core_tests PROPERTIES SKIP_RETURN_CODE 4)
