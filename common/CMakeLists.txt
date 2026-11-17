# Header-only definitions shared by the simulation and the analysis: units, the
# on-disk event schema and the data-directory policy. This is the only contract
# between sim/ and analysis/.
configure_file(include/riptide/common/version.hpp.in
    "${CMAKE_CURRENT_BINARY_DIR}/include/riptide/common/version.hpp" @ONLY)

# git_info.hpp is regenerated at every build (see cmake/RiptideGitDescribe.cmake);
# it is also generated once here so the header exists before the first build.
set(riptide_git_info_header "${CMAKE_CURRENT_BINARY_DIR}/include/riptide/common/git_info.hpp")
set(riptide_git_info_args
    "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}"
    "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DTEMPLATE=${CMAKE_CURRENT_SOURCE_DIR}/include/riptide/common/git_info.hpp.in"
    "-DOUTPUT=${riptide_git_info_header}")
set(riptide_git_info_script "${PROJECT_SOURCE_DIR}/cmake/RiptideGitDescribe.cmake")
execute_process(COMMAND "${CMAKE_COMMAND}" ${riptide_git_info_args} -P "${riptide_git_info_script}"
    COMMAND_ERROR_IS_FATAL ANY)
add_custom_target(riptide_git_info ALL
    COMMAND "${CMAKE_COMMAND}" ${riptide_git_info_args} -P "${riptide_git_info_script}"
    BYPRODUCTS "${riptide_git_info_header}"
    COMMENT "Recording git describe"
    VERBATIM)

add_library(riptide_common INTERFACE)
add_library(riptide::common ALIAS riptide_common)
target_include_directories(riptide_common INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")
target_compile_features(riptide_common INTERFACE cxx_std_17)
# Dependencies of an INTERFACE library are followed by its consumers, so every
# target using riptide::common is built after the header is refreshed.
add_dependencies(riptide_common riptide_git_info)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
