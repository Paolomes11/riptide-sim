# Vendored header-only dependencies (see each VERSION file for source and
# sha256). Include directories are SYSTEM so that the project warnings and
# -Werror do not fire inside third-party headers.

add_library(riptide_ext_json INTERFACE)
add_library(riptide::ext_json ALIAS riptide_ext_json)
target_include_directories(riptide_ext_json SYSTEM INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/nlohmann_json/include")

add_library(riptide_ext_lyra INTERFACE)
add_library(riptide::ext_lyra ALIAS riptide_ext_lyra)
target_include_directories(riptide_ext_lyra SYSTEM INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/lyra/include")
