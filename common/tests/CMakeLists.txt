add_executable(common_tests
    data_paths_test.cpp
    event_sample_test.cpp
    event_schema_test.cpp
    hp_kinematics_test.cpp
    np_kinematics_test.cpp
    random_test.cpp
    units_test.cpp)
target_link_libraries(common_tests PRIVATE
    riptide::common riptide::options riptide::warnings Catch2::Catch2WithMain)
# SKIP_RETURN_CODE 4: Catch2 exits with 4 when every selected test was SKIP()ped.
catch_discover_tests(common_tests PROPERTIES SKIP_RETURN_CODE 4)
