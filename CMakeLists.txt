cmake_minimum_required(VERSION 3.10)
project(MyProject C)

# 1. Compiler Flags
set(CMAKE_C_STANDARD 17)
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
    -Walloca
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

include(CTest)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
# Define the internal libraries group
set(INTERNAL_LIBS
    mylib
    mysound
    mysync
    mysignal
    myfft
    mymath
    myarr
)

set(SRC ${PROJECT_ROOT}/src)
# 3. Create Internal Static Libraries
# Equivalent to $(LIB)/libmylib.a
add_library(mylib STATIC ${SRC}/mylib.c ${SRC}/crosscorint.c ${SRC}/parseargs.c)
target_link_libraries(mylib PUBLIC myarr )
add_library(mysound STATIC ${SRC}/mysound.c)
target_link_libraries(mysound PUBLIC myarr )
add_library(mysync STATIC ${SRC}/mysync.c)
target_link_libraries(mysync PUBLIC myarr )
add_library(mysignal STATIC ${SRC}/mysignal.c)
add_library(myfft STATIC ${SRC}/myfft.c)
target_link_libraries(myfft PUBLIC mysync myarr )
add_library(mymath STATIC ${SRC}/mymath.c)
target_link_libraries(mymath PUBLIC myarr )
add_library(myarr STATIC ${SRC}/myarr.c)


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
    target_link_libraries(${BIN} 
        ${INTERNAL_LIBS} 
        ${ASOUND_LIB} 
        ${FFTW3_LIB} 
        ${MATH_LIB}
    )
endforeach()


# 6. Install (Moves binaries to a bin folder)
install(TARGETS ${BINARIES} DESTINATION bin)

