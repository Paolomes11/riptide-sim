add_executable(sim_config_tests
    cli_options_test.cpp
    event_seeds_test.cpp
    run_config_test.cpp
    run_record_test.cpp)
target_link_libraries(sim_config_tests PRIVATE
    riptide::sim_config riptide::options riptide::warnings Catch2::Catch2WithMain)
# The shipped config files are parsed by the tests, so they cannot rot unnoticed.
target_compile_definitions(sim_config_tests PRIVATE
    RIPTIDE_CONFIG_DIR="${PROJECT_SOURCE_DIR}/config")
catch_discover_tests(sim_config_tests PROPERTIES SKIP_RETURN_CODE 4)
