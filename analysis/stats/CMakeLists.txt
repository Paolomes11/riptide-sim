add_library(riptide_stats STATIC
    src/case_a.cpp
    src/combine.cpp
    src/grid.cpp
    src/numerics.cpp
    src/priors.cpp)
add_library(riptide::stats ALIAS riptide_stats)
target_include_directories(riptide_stats PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(riptide_stats
    PUBLIC riptide::common
    PRIVATE riptide::options riptide::warnings)

if(RIPTIDE_BUILD_TESTS)
    add_subdirectory(tests)
endif()
