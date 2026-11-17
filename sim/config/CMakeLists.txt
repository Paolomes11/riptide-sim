# Run configuration (JSON -> RunConfig). Geant4-free on purpose: it is added
# from the top level even when RIPTIDE_BUILD_SIM is OFF, so the parser and its
# validation are tested in the CI build too.
add_library(riptide_sim_config STATIC
    src/run_config.cpp)
add_library(riptide::sim_config ALIAS riptide_sim_config)
target_include_directories(riptide_sim_config PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_sim_config
    PUBLIC riptide::common riptide::ext_json
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
