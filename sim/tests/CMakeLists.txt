add_executable(sim_tests
    units_consistency_test.cpp)
target_link_libraries(sim_tests PRIVATE
    riptide::common riptide::options riptide::warnings Catch2::Catch2WithMain ${Geant4_LIBRARIES})
catch_discover_tests(sim_tests PROPERTIES SKIP_RETURN_CODE 4)
