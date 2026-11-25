add_library(riptide_stats STATIC
    src/calibration.cpp
    src/acceptance.cpp
    src/case_a.cpp
    src/case_b.cpp
    src/case_b_posterior.cpp
    src/pseudo_experiment_b.cpp
    src/combine.cpp
    src/direction_grid.cpp
    src/vmf_kernel.cpp
    src/grid.cpp
    src/numerics.cpp
    src/priors.cpp
    src/pseudo_experiment.cpp src/split.cpp
    src/validate.cpp)
add_library(riptide::stats ALIAS riptide_stats)
target_include_directories(riptide_stats PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_stats
    PUBLIC riptide::common
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
