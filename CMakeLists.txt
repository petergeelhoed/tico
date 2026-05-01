cmake_minimum_required(VERSION 3.10)
project(MyProject C)

# 1. Compiler Flags
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)



# =========================================
# Base: portable across Clang & GCC (C/C++)
# =========================================
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Werror

    # Your existing set
    -Wshadow
    -Wpointer-arith
    -Wcast-qual
    -Wcast-align
    -Wmissing-prototypes          # C only
    -Wstrict-prototypes           # C only
    -Wconversion
    -Wsign-conversion
    -Wformat=2
    -Wundef
    -Wwrite-strings
    -funsigned-char
    -Wswitch-enum
    -Wswitch-default
    -Wimplicit-fallthrough
    -Wmissing-field-initializers
    -Wnull-dereference
    -Wdouble-promotion
    -Wold-style-definition        # C only
    -Wmissing-declarations
    -Wredundant-decls
    -Wbad-function-cast           # C only (GCC/Clang both accept)
    -Wnested-externs              # C only
    -Wformat-security
    -Wdate-time
    -Wsign-compare
    -Warray-bounds                # (no "=2" in Clang)
    -Wvla
#   -Walloca
    -Wfloat-equal                 # NaN/precision foot-guns
    -Wtautological-compare        # obvious always true/false
    -Woverlength-strings          # extremely long string literals
)

# ======================================================
# Clang-specific / Clang-strong warnings (Clang 15+)
# ======================================================
if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  add_compile_options(
    # You had these:
    -Wdocumentation
    -Wunreachable-code
    -Wformat-nonliteral

    # More Clang goodness:
    -Wcomma                               # comma used where logical ops intended
    -Wbitwise-instead-of-logical          # a | b used instead of a || b
    -Wlogical-not-parentheses             # ! applied without parentheses
    -Wimplicit-int                        # C: missing type defaults to int (now error-prone in Clang 15)
    -Wimplicit-function-declaration       # C: implicit decls (tightened in Clang 15)
    -Wdeprecated-non-prototype            # C: decls without prototype (C2x semantics)
    -Wreserved-identifier                 # identifiers starting with _ or __
    -Wnewline-eof                         # final newline hygiene
    -Warray-parameter                     # array-to-pointer decay surprises in prototypes
    -Warray-bounds-pointer-arithmetic     # pointer arithmetic going past arrays
    -Wassign-enum                         # assigning ints to enums etc.
    -Wswitch-bool                         # switch on boolean suspicious
    -Wshadow-uncaptured-local     # C++ lambdas; harmless on C
    -Wconditional-uninitialized   # uninit via conditionals (Clang strong)
  )
endif()


# ==========================
# GCC-specific (your list)
# ==========================
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  add_compile_options(
    -Wrestrict
    -Wstringop-overflow=4
    -Wstringop-truncation
    -Wformat-overflow=2
    -Wformat-truncation=2
    -Wshift-overflow=2
    -Wvla-larger-than=8192
    -Walloca-larger-than=1024
    -Wstrict-overflow=3
    -Warray-bounds=2
    -Wlogical-op                        # GCC’s version of logical-op pitfalls
    -Wduplicated-cond                   # duplicated conditions in if/else if
    -Wduplicated-branches               # duplicated branches
  )
endif()




# Prefer portable thread discovery; this also adds the correct link flags.
find_package(Threads REQUIRED)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Sanitizer Options (Uncomment one to enable)
# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address -static-libasan)

# 2. External Library Dependencies
find_library(ASOUND_LIB asound)
find_library(FFTW3_LIB fftw3)
find_library(MATH_LIB m)

find_program(CLANG_FORMAT_EXE
  NAMES clang-format clang-format-19 clang-format-18 clang-format-17 clang-format-16 clang-format-15)
find_program(CLANG_TIDY_EXE
  NAMES clang-tidy clang-tidy-19 clang-tidy-18 clang-tidy-17 clang-tidy-16 clang-tidy-15)

file(GLOB_RECURSE FORMAT_SOURCES CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/src/*.c
  ${PROJECT_SOURCE_DIR}/inc/*.h
  ${PROJECT_SOURCE_DIR}/tests/*.c)

file(GLOB_RECURSE TIDY_SOURCES CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/src/*.c
  ${PROJECT_SOURCE_DIR}/inc/*.h)

if (CLANG_FORMAT_EXE)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXE} -i ${FORMAT_SOURCES}
    COMMENT "Formatting C sources with clang-format")

  add_custom_target(format-check
    COMMAND ${CLANG_FORMAT_EXE} --dry-run --Werror ${FORMAT_SOURCES}
    COMMENT "Checking C sources formatting with clang-format")
else()
  add_custom_target(format
    COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
    COMMAND ${CMAKE_COMMAND} -E false
    COMMENT "clang-format is required for format target")

  add_custom_target(format-check
    COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
    COMMAND ${CMAKE_COMMAND} -E false
    COMMENT "clang-format is required for format-check target")
endif()

if (CLANG_TIDY_EXE)
  add_custom_target(clang-tidy
  COMMAND ${CLANG_TIDY_EXE}
    -p ${CMAKE_BINARY_DIR}
    --extra-arg=-Wno-unknown-warning-option
    --extra-arg=-Wno-error=unknown-warning-option
    ${TIDY_SOURCES}
    COMMENT "Running clang-tidy over C sources")
else()
  add_custom_target(clang-tidy
    COMMAND ${CMAKE_COMMAND} -E echo "clang-tidy not found"
    COMMAND ${CMAKE_COMMAND} -E false
    COMMENT "clang-tidy is required for clang-tidy target")
endif()

# Runs a GCC debug build from the source root so quality can include GCC compilation.
add_custom_target(quality-gcc-build
  COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR} ${CMAKE_COMMAND} --preset debug-gcc
  COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR} ${CMAKE_COMMAND} --build --preset debug-gcc
  COMMENT "Configuring and building GCC debug preset")

include(CTest)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
# Define the internal libraries group
set(INTERNAL_LIBS
  analysis
  capture_helpers
  cli
    mysound
    mysync
    mysignal
    myfft
    mymath
    myarr
)

set(SRC ${PROJECT_ROOT}/src)
# 3. Create Internal Static Libraries

add_library(myarr STATIC ${SRC}/myarr.c)

add_library(mysync STATIC ${SRC}/mysync.c)
target_link_libraries(mysync PUBLIC myarr)

add_library(mysignal STATIC ${SRC}/mysignal.c)

add_library(mymath STATIC ${SRC}/mymath.c)
target_link_libraries(mymath PUBLIC myarr)

add_library(myfft STATIC ${SRC}/myfft.c)
target_link_libraries(myfft PUBLIC mysync myarr)

add_library(mysound STATIC ${SRC}/mysound.c)
target_link_libraries(mysound PUBLIC myarr)

add_library(cli STATIC ${SRC}/parseargs.c)

add_library(analysis STATIC ${SRC}/analysis.c ${SRC}/crosscorint.c)
target_link_libraries(analysis PUBLIC myarr mymath mysync myfft cli)

add_library(capture_helpers STATIC ${SRC}/capture_helpers.c ${SRC}/resources.c)
target_link_libraries(capture_helpers
    PUBLIC myarr mysound analysis myfft mymath mysync cli)


foreach(LIB ${INTERNAL_LIBS})
    target_include_directories(${LIB} PUBLIC ${PROJECT_ROOT}/inc)
    target_link_libraries(${LIB} PUBLIC Threads::Threads)
endforeach()


# 5. Define Test Targets
set(BINARIES
    capture
    1dcor
    fft
    linregw
    ifft
    linreg
    record
)
add_subdirectory(tests)

foreach(BIN ${BINARIES})
    add_executable(${BIN} src/${BIN}.c)
endforeach()

target_link_libraries(capture
  capture_helpers
  mysignal
  mysound
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(1dcor
  myfft
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(fft
  mymath
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(linregw
  mymath
  cli
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(ifft
  cli
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(linreg
  cli
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})

target_link_libraries(record
  mysound
  mysync
  cli
  myarr
  ${ASOUND_LIB}
  ${FFTW3_LIB}
  ${MATH_LIB})


# 6. Install (Moves binaries to a bin folder)
install(TARGETS ${BINARIES} DESTINATION bin)

