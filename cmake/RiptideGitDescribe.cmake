# Script mode (cmake -P): writes the output of `git describe` into a header, so
# every run can record the exact source it was built from.
#
# Inputs (-D): GIT_EXECUTABLE (may be empty), SOURCE_DIR, TEMPLATE, OUTPUT.
#
# It runs at every build, not only at configure time: a value captured by
# `cmake --preset` would go stale after the next commit. configure_file leaves
# OUTPUT untouched when its content is unchanged, so an unchanged commit does
# not trigger any recompilation.
#
# The value describes the source tree at build time; "-dirty" marks
# uncommitted changes. Without git or outside a repository it is "unknown".
set(RIPTIDE_GIT_DESCRIBE "unknown")
if(GIT_EXECUTABLE)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --tags --always --dirty
        WORKING_DIRECTORY "${SOURCE_DIR}"
        RESULT_VARIABLE git_result
        OUTPUT_VARIABLE git_output
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(git_result EQUAL 0 AND NOT git_output STREQUAL "")
        set(RIPTIDE_GIT_DESCRIBE "${git_output}")
    endif()
endif()
configure_file("${TEMPLATE}" "${OUTPUT}" @ONLY)
