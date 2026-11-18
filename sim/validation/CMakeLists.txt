# Physics checks of riptide_sim against evaluated nuclear data: ENDF reading,
# expected distributions and the statistics of the comparisons. Geant4-free
# on purpose, like sim/config and sim/truth: it is added from the top level
# and tested in every preset; the Geant4 app that reads the ntuples links it.
add_library(riptide_sim_validation STATIC
    src/endf_format.cpp src/endf_reader.cpp src/check_stats.cpp
    src/expected_distributions.cpp src/sha256.cpp src/beam_geometry.cpp)
add_library(riptide::sim_validation ALIAS riptide_sim_validation)
target_include_directories(riptide_sim_validation PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_sim_validation
    PUBLIC riptide::common riptide::sim_config
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
