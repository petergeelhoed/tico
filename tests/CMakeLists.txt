
# tests/CMakeLists.txt

set(TEST_TARGETS 
    testfft
    testsyncwrite
    testsyncappend
    testgetdoubles
    testmysync testfilter testremove50 testlinreg 
    testinvert testmulmat testmatlinreg testmatlinregquad 
    testcalculateall testcrosscorint testtranspone
    testshiftbufferdata testrotatederivativewindow testupdatetotalshift
    testcheckargs testgetbeaterror
    testgetinput
    testmicamplification
)
foreach(TEST ${TEST_TARGETS})
    add_executable(${TEST} ${TEST}.c ${CMAKE_SOURCE_DIR}/src/columns.c)
    target_include_directories(${TEST} PRIVATE ${CMAKE_SOURCE_DIR}/inc)
endforeach()

target_link_libraries(testfft myfft ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testsyncwrite myarr mysync ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testsyncappend myarr mysync ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testgetdoubles cli ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testmysync myarr mysync ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testfilter myfft ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testremove50 myfft ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testlinreg mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testinvert mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testmulmat mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testmatlinreg mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testmatlinregquad mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testcalculateall analysis ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testcrosscorint analysis ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testtranspone mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testgetinput cli mymath ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testmicamplification mysound capture_helpers printing ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testshiftbufferdata ${INTERNAL_LIBS} capture_helpers printing ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testrotatederivativewindow ${INTERNAL_LIBS} capture_helpers printing ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testupdatetotalshift ${INTERNAL_LIBS} capture_helpers printing ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testcheckargs cli ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})
target_link_libraries(testgetbeaterror analysis ${ASOUND_LIB} ${FFTW3_LIB} ${MATH_LIB})

target_compile_definitions(testcalculateall
    PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set( TESTS
        testcalculateall
        testcrosscorint
        testinvert
        testshiftbufferdata
        testrotatederivativewindow
        testupdatetotalshift
        testcheckargs
        testgetbeaterror
        )
foreach(TST ${TESTS})

if(TST STREQUAL "testcalculateall")
    add_test(NAME ${TST} COMMAND $<TARGET_FILE:${TST}>)
elseif(TST STREQUAL "testshiftbufferdata" OR
       TST STREQUAL "testrotatederivativewindow" OR
    TST STREQUAL "testupdatetotalshift" OR
    TST STREQUAL "testcheckargs" OR
    TST STREQUAL "testgetbeaterror")
    add_test(NAME ${TST} COMMAND $<TARGET_FILE:${TST}>)
else()
    add_test(
      NAME ${TST}
      COMMAND /bin/sh -c
              "set -e; \
               \"$<TARGET_FILE:${TST}>\" > \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\"; \
               \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_${TST}.txt\""
    )
endif()

endforeach()
