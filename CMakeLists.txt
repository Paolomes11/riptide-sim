cmake_minimum_required(VERSION 3.25)

project(riptide_sim
    VERSION 0.5.0
    DESCRIPTION "Geant4 simulation and C++ statistical reconstruction of incident neutrons from n-p scattering"
    LANGUAGES CXX)

# ROOT 6.36 on the target machine is built with C++17: every component that may
# share headers or objects with ROOT must use the same standard to avoid ABI and
# dictionary mismatches.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(RIPTIDE_BUILD_SIM "Build the Geant4 application (requires Geant4)" ON)
option(RIPTIDE_BUILD_ROOT_IO "Build the ROOT reader of the analysis (requires ROOT 6.36)" ON)
option(RIPTIDE_BUILD_TESTS "Build the Catch2 unit tests" ON)
option(RIPTIDE_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RIPTIDE_ENABLE_SANITIZERS "Build with AddressSanitizer and UBSan" OFF)
option(RIPTIDE_ENABLE_NATIVE "Optimize for the build machine (-march=native); binaries are not portable" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(RiptideCompilerOptions)

# Optional: without git the build still works and records "unknown" as the
# source revision (cmake/RiptideGitDescribe.cmake).
find_package(Git QUIET)

if(RIPTIDE_BUILD_TESTS)
    enable_testing()
    find_package(Catch2 3 REQUIRED)
    include(Catch)
endif()

# Dependency direction (only downwards):
#   common  <- sim
#   common  <- analysis/stats   (pure C++: no ROOT, no Geant4)
#   common  <- analysis/io/core <- analysis/io/root (the only ROOT user)
# sim/config, sim/truth and sim/validation have no Geant4 dependency and are
# built (and tested) in every preset.
add_subdirectory(external)
add_subdirectory(common)
add_subdirectory(sim/config)
add_subdirectory(sim/truth)
add_subdirectory(sim/validation)
add_subdirectory(analysis)
if(RIPTIDE_BUILD_SIM)
    add_subdirectory(sim)
endif()
