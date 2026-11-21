# Run configuration (JSON -> RunConfig) and command line of riptide_sim.
# Geant4-free on purpose: it is added from the top level even when
# RIPTIDE_BUILD_SIM is OFF, so parsing and validation are tested in CI too.
add_library(riptide_sim_config STATIC
    src/cli_options.cpp
    src/run_config.cpp
    src/run_record.cpp
    src/source_sampling.cpp)
add_library(riptide::sim_config ALIAS riptide_sim_config)
target_include_directories(riptide_sim_config PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_sim_config
    PUBLIC riptide::common riptide::ext_json
    PRIVATE riptide::ext_lyra riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
