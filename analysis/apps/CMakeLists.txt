# The analysis applications. Their command lines are a pure library, built and
# tested in every preset; the executables read the event files and so are
# built only with the ROOT reader.
add_library(riptide_analysis_cli STATIC src/prepare_cli.cpp)
add_library(riptide::analysis_cli ALIAS riptide_analysis_cli)
target_include_directories(riptide_analysis_cli PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_analysis_cli
    PUBLIC riptide::io_core
    PRIVATE riptide::ext_lyra riptide::options riptide::warnings)

# The pure part of riptide_case_a_geant4: no ROOT, built and tested everywhere.
add_library(riptide_case_a_geant4_core STATIC src/case_a_geant4.cpp)
add_library(riptide::case_a_geant4_core ALIAS riptide_case_a_geant4_core)
target_include_directories(riptide_case_a_geant4_core PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_case_a_geant4_core
    PUBLIC riptide::analysis_cli riptide::stats riptide::io_core riptide::common
    PRIVATE riptide::options riptide::warnings)

# The statistics applications need neither ROOT nor Geant4: built everywhere.
add_executable(riptide_stats_bench src/riptide_stats_bench.cpp)
target_link_libraries(riptide_stats_bench PRIVATE
    riptide::stats riptide::common riptide::ext_lyra riptide::options riptide::warnings)
add_executable(riptide_stats_calibrate src/riptide_stats_calibrate.cpp)
target_link_libraries(riptide_stats_calibrate PRIVATE
    riptide::case_a_geant4_core riptide::stats riptide::io_core riptide::common riptide::ext_json riptide::ext_lyra
    riptide::options riptide::warnings)

if(RIPTIDE_BUILD_ROOT_IO)
    add_executable(riptide_analysis_prepare src/riptide_analysis_prepare.cpp)
    target_link_libraries(riptide_analysis_prepare PRIVATE
        riptide::analysis_cli riptide::io_core riptide::io_root riptide::common
        riptide::options riptide::warnings)
    add_executable(riptide_case_a_geant4 src/riptide_case_a_geant4.cpp)
    target_link_libraries(riptide_case_a_geant4 PRIVATE
        riptide::case_a_geant4_core riptide::io_core riptide::io_root riptide::stats
        riptide::common riptide::ext_json riptide::ext_lyra riptide::options riptide::warnings)
    add_executable(riptide_acceptance src/riptide_acceptance.cpp)
    target_link_libraries(riptide_acceptance PRIVATE
        riptide::analysis_cli riptide::io_core riptide::io_root riptide::common
        riptide::ext_json riptide::ext_lyra riptide::options riptide::warnings)
endif()

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
