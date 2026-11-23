cmake_minimum_required(VERSION 3.25)
project(dop14-area LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

enable_testing()

# CornerCut<f32> must round as written: no FMA contraction
add_library(dop14_flags INTERFACE)
target_compile_options(dop14_flags INTERFACE
    -ffp-contract=off
    -fno-exceptions
    -Wall -Wextra -Wpedantic -Wshadow
)

# system cddlib, GMP build (Arch: cddlib, Debian/Ubuntu: libcdd-dev); headers
# in cddlib/ upstream, cdd/ on some distributions
add_library(dop14_cddgmp INTERFACE)
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
    pkg_check_modules(CDDGMP QUIET IMPORTED_TARGET cddgmp)
endif()
find_path(CDD_INCLUDE_DIR setoper.h
    HINTS ${CDDGMP_INCLUDEDIR} ${CDDGMP_INCLUDE_DIRS}
    PATH_SUFFIXES cddlib cdd
    REQUIRED
)
if (CDDGMP_FOUND)
    target_link_libraries(dop14_cddgmp INTERFACE PkgConfig::CDDGMP)
else()
    find_library(CDDGMP_LIBRARY cddgmp REQUIRED)
    find_library(GMP_LIBRARY gmp REQUIRED)
    target_link_libraries(dop14_cddgmp INTERFACE ${CDDGMP_LIBRARY} ${GMP_LIBRARY})
    target_compile_definitions(dop14_cddgmp INTERFACE GMPRATIONAL)
endif()
target_include_directories(dop14_cddgmp SYSTEM INTERFACE ${CDD_INCLUDE_DIR})
message(STATUS "cddlib (GMP): headers in ${CDD_INCLUDE_DIR}, pkg-config cddgmp ${CDDGMP_VERSION}")

add_library(dop14_core STATIC
    src/cases.cpp
    src/ref_cdd.cpp
    src/ref_brute.cpp
)
target_include_directories(dop14_core PUBLIC src)
target_link_libraries(dop14_core PUBLIC dop14_flags PRIVATE dop14_cddgmp)

# timing.cpp calls cddlib directly
add_executable(dop14-area src/main.cpp src/timing.cpp)
target_link_libraries(dop14-area PRIVATE dop14_core dop14_cddgmp)

add_subdirectory(tests)
