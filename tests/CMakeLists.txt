
# tests/CMakeLists.txt

set(TEST_TARGETS 
    testfft
    testgetdoubles
    testmysync testfilter testremove50 testlinreg 
    testinvert testmulmat testmatlinreg testmatlinregquad 
    testcalculateall testcrosscorint testtranspone
    testgetinput
)
foreach(TEST ${TEST_TARGETS})
    add_executable(${TEST} ${TEST}.c)
    target_link_libraries(${TEST} 
        ${INTERNAL_LIBS} 
        ${ASOUND_LIB} 
        ${FFTW3_LIB} 
        ${MATH_LIB}
    )
endforeach()

set( TESTS
        testcrosscorint
        testinvert 
        )
foreach(TST ${TESTS})

add_test(
  NAME ${TST}
  COMMAND /bin/sh -c
          "set -e; \
           \"$<TARGET_FILE:${TST}>\" > \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\"; \
           \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_${TST}.txt\""
)

endforeach()
