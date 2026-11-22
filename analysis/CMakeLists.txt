# analysis/stats and analysis/io/core are pure C++ (no ROOT, no Geant4) so they
# build and test anywhere, including CI. analysis/io/root is the only part that
# needs ROOT; the applications come later.
add_subdirectory(stats)
add_subdirectory(io/core)
if(RIPTIDE_BUILD_ROOT_IO)
    add_subdirectory(io/root)
endif()
