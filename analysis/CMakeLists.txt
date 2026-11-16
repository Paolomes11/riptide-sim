# analysis/stats is pure C++ (no ROOT, no Geant4) so it builds and tests anywhere,
# including CI. The ROOT reader (analysis/io) and the applications come later.
add_subdirectory(stats)
