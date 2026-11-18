# Reader of riptide_sim runs (run.json plus the per-thread ntuples), shared by
# the integration test and the phase-1b check app. It reads through the Geant4
# ROOT reader, so it needs Geant4 and is added from sim/, not from the top
# level; its public header stays Geant4-free.
add_library(riptide_sim_reader STATIC
    src/run_reader.cpp)
add_library(riptide::sim_reader ALIAS riptide_sim_reader)
target_include_directories(riptide_sim_reader PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_sim_reader
    PUBLIC riptide::common riptide::sim_config riptide::sim_truth
    PRIVATE riptide::options riptide::warnings ${Geant4_LIBRARIES})

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
