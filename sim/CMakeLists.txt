# Geant4 application. Geant4 is found through CMAKE_PREFIX_PATH or PATH (the
# directory holding geant4-config); no install path is hard-coded here.
find_package(Geant4 11.4 REQUIRED)

add_executable(riptide_sim
    apps/riptide_sim.cpp
    src/action_initialization.cpp
    src/detector_construction.cpp
    src/event_action.cpp
    src/primary_generator_action.cpp
    src/run_action.cpp
    src/stepping_action.cpp)
target_include_directories(riptide_sim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
# Imported Geant4 targets carry their include directories as SYSTEM, so the
# project warnings do not fire inside Geant4 headers.
target_link_libraries(riptide_sim PRIVATE
    riptide::common riptide::sim_config riptide::sim_truth riptide::options riptide::warnings ${Geant4_LIBRARIES})

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
