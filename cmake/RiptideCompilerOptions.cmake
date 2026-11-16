# Project-wide compiler settings, applied per target through two INTERFACE
# targets so that third-party code (Geant4, ROOT, Catch2) is never affected.

add_library(riptide_warnings INTERFACE)
add_library(riptide::warnings ALIAS riptide_warnings)
target_compile_options(riptide_warnings INTERFACE
    -Wall -Wextra -Wpedantic
    -Wshadow -Wconversion -Wsign-conversion -Wdouble-promotion
    -Wold-style-cast -Wcast-qual -Wcast-align
    -Wnon-virtual-dtor -Woverloaded-virtual
    -Wnull-dereference -Wimplicit-fallthrough -Wformat=2 -Wundef
    -Wextra-semi -Wunused
    # exact float comparison is almost always a bug in numerical code
    -Werror=float-equal
    $<$<BOOL:${RIPTIDE_WARNINGS_AS_ERRORS}>:-Werror>)

add_library(riptide_options INTERFACE)
add_library(riptide::options ALIAS riptide_options)
# -ffast-math is deliberately absent: it lets the compiler assume no NaN/Inf and
# reorder sums, which breaks log-likelihoods that are legitimately -inf and the
# reproducibility of the statistical results.
target_compile_options(riptide_options INTERFACE
    $<$<BOOL:${RIPTIDE_ENABLE_NATIVE}>:-march=native>)

if(RIPTIDE_ENABLE_SANITIZERS)
    target_compile_options(riptide_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_link_options(riptide_options INTERFACE -fsanitize=address,undefined)
endif()
