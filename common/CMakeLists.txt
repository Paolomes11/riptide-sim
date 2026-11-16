# Header-only definitions shared by the simulation and the analysis: units, the
# on-disk event schema and the data-directory policy. This is the only contract
# between sim/ and analysis/.
configure_file(include/riptide/common/version.hpp.in
    "${CMAKE_CURRENT_BINARY_DIR}/include/riptide/common/version.hpp" @ONLY)

add_library(riptide_common INTERFACE)
add_library(riptide::common ALIAS riptide_common)
target_include_directories(riptide_common INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")
target_compile_features(riptide_common INTERFACE cxx_std_17)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
