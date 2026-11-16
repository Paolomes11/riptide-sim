# Geant4 application. Geant4 is found through CMAKE_PREFIX_PATH or PATH (the
# directory holding geant4-config); no install path is hard-coded here.
find_package(Geant4 11.4 REQUIRED)

add_executable(riptide_sim apps/riptide_sim.cpp)
# Imported Geant4 targets carry their include directories as SYSTEM, so the
# project warnings do not fire inside Geant4 headers.
target_link_libraries(riptide_sim PRIVATE
    riptide::common riptide::options riptide::warnings ${Geant4_LIBRARIES})

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
