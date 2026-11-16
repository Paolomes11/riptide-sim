// Entry point of the Geant4 simulation. Phase 0 skeleton: it only reports the
// versions it was built with, to check the toolchain end to end.

#include <G4Version.hh>

#include <iostream>
#include <string_view>

#include "riptide/common/event_schema.hpp"
#include "riptide/common/version.hpp"

namespace {

void printVersion() {
    std::cout << "riptide_sim " << riptide::kProjectVersion << '\n'
              << "Geant4 " << G4VERSION_NUMBER / 100 << '.' << (G4VERSION_NUMBER / 10) % 10 << '.'
              << G4VERSION_NUMBER % 10 << ' ' << G4Version << ' ' << G4Date << '\n'
              << "event schema " << riptide::schema::kSchemaVersion << '\n';
}

void printUsage() {
    std::cout << "usage: riptide_sim --version\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view{argv[1]} == "--version") {
        printVersion();
        return 0;
    }
    printUsage();
    return argc == 1 ? 0 : 2;
}
