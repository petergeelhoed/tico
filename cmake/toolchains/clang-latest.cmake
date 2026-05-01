# Pick the newest available clang-X binary first, then fallback to plain clang.
set(_clang_candidates
    clang-24
    clang-23
    clang-22
    clang-21
    clang-20
    clang-19
    clang-18
    clang-17
    clang-16
    clang-15
    clang-14
    clang)

find_program(CLANG_LATEST_COMPILER NAMES ${_clang_candidates})

if(NOT CLANG_LATEST_COMPILER)
    message(FATAL_ERROR "No clang compiler found. Install clang or add it to PATH.")
endif()

set(CMAKE_C_COMPILER "${CLANG_LATEST_COMPILER}" CACHE FILEPATH "Detected latest clang compiler" FORCE)