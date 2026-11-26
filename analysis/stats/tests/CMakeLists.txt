add_executable(stats_tests
    acceptance_test.cpp
    calibration_test.cpp
    calibration_b_test.cpp
    case_a_test.cpp
    case_b_test.cpp
    case_b_posterior_test.cpp
    case_c_test.cpp
    case_c_posterior_test.cpp
    pseudo_experiment_b_test.cpp
    pseudo_experiment_c_test.cpp
    resample_c_test.cpp
    calibration_c_test.cpp
    direction_grid_test.cpp
    vmf_kernel_test.cpp
    grid_test.cpp
    log_sum_exp_test.cpp
    numerics_test.cpp
    priors_test.cpp
    pseudo_experiment_test.cpp
    split_test.cpp
    validate_test.cpp)
target_link_libraries(stats_tests PRIVATE
    riptide::stats riptide::options riptide::warnings Catch2::Catch2WithMain)
catch_discover_tests(stats_tests PROPERTIES SKIP_RETURN_CODE 4)
