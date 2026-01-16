cmake_minimum_required(VERSION 3.10)
project(MyProject C)

# 1. Compiler Flags
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

add_compile_options(-Wall -Wpedantic -Wextra -Wsign-compare -Werror -g)

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
add_library(mysound STATIC ${SRC}/mysound.c)
add_library(mysync STATIC ${SRC}/mysync.c)
target_link_libraries(mysync PUBLIC myarr )
add_library(mysignal STATIC ${SRC}/mysignal.c)
add_library(myfft STATIC ${SRC}/myfft.c)
target_link_libraries(myfft PUBLIC mysync )
add_library(mymath STATIC ${SRC}/mymath.c)
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

