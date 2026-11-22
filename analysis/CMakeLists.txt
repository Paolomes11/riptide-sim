# analysis/stats and analysis/io/core are pure C++ (no ROOT, no Geant4) so they
# build and test anywhere, including CI. The ROOT reader and the applications
# come later.
add_subdirectory(stats)
add_subdirectory(io/core)
