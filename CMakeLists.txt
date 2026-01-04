cmake_minimum_required(VERSION 3.10)
project(MyProject C)

# 1. Compiler Flags
set(CMAKE_C_STANDARD 99)
add_compile_options(-Wall -Wpedantic -Wextra -Wsign-compare -Werror -pthread -g)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Sanitizer Options (Uncomment one to enable)
# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address -static-libasan)

# 2. External Library Dependencies
find_library(ASOUND_LIB asound)
find_library(FFTW3_LIB fftw3)
find_library(MATH_LIB m)

# 3. Create Internal Static Libraries
# Equivalent to $(LIB)/libmylib.a
add_library(mylib STATIC mylib.c crosscorint.c parseargs.c)
add_library(mysound STATIC mysound.c)
add_library(mysync STATIC mysync.c)
add_library(mysignal STATIC mysignal.c)
add_library(myfft STATIC myfft.c)
target_link_libraries(myfft PRIVATE mysync)
add_library(mymath STATIC mymath.c)

# 5. Define Test Targets
set(TEST_TARGETS 
 testfft
 capture 1dcor
    testmysync testfilter testremove50 testlinreg 
    testinvert testmulmat testmatlinreg testmatlinregquad 
    testcalculateall testcrosscorint testtranspone
)
# Define the internal libraries group
set(INTERNAL_LIBS 
    mylib 
    mysound 
    mysync 
    mysignal 
    myfft 
    mymath
)

# Use -Wl,--start-group and -Wl,--end-group to handle circular dependencies
foreach(TGT ${TARGETS})
    add_executable(${TGT} ${TGT}.c)
    target_link_libraries(${TGT} 
        ${INTERNAL_LIBS} 
        ${ASOUND_LIB} 
        ${FFTW3_LIB} 
        ${MATH_LIB}
    )
endforeach()

foreach(TEST ${TEST_TARGETS})
    add_executable(${TEST} ${TEST}.c myarr.c)
    target_link_libraries(${TEST} 
        ${INTERNAL_LIBS} 
        ${ASOUND_LIB} 
        ${FFTW3_LIB} 
        ${MATH_LIB}
    )
endforeach()


# 6. Install (Moves binaries to a bin folder)
install(TARGETS ${TARGETS} DESTINATION bin)

