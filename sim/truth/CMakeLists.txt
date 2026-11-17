# Monte Carlo truth bookkeeping of riptide_sim. Geant4-free on purpose, like
# sim/config: it is added from the top level and tested in every preset.
add_library(riptide_sim_truth STATIC
    src/truth_record.cpp)
add_library(riptide::sim_truth ALIAS riptide_sim_truth)
target_include_directories(riptide_sim_truth PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_sim_truth
    PUBLIC riptide::common riptide::sim_config
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
