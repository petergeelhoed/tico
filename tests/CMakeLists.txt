
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
    testerf
)
foreach(TEST ${TEST_TARGETS})
    add_executable(${TEST} ${TEST}.c)
    target_include_directories(${TEST} PRIVATE ${CMAKE_SOURCE_DIR}/inc)
endforeach()

target_link_libraries(testfft myfft )
target_link_libraries(testsyncwrite mysync)
target_link_libraries(testsyncappend mysync)
target_link_libraries(testgetdoubles cli)
target_link_libraries(testmysync mysync)
target_link_libraries(testfilter myfft)
target_link_libraries(testremove50 myfft)
target_link_libraries(testlinreg mymath )
target_link_libraries(testinvert mymath )
target_link_libraries(testmulmat mymath )
target_link_libraries(testmatlinreg mymath )
target_link_libraries(testmatlinregquad mymath)
target_link_libraries(testcalculateall analysis)
target_link_libraries(testcrosscorint analysis )
target_link_libraries(testtranspone mymath)
target_link_libraries(testgetinput cli m)
target_link_libraries(testmicamplification mysound )
target_link_libraries(testshiftbufferdata capture_helpers)
target_link_libraries(testrotatederivativewindow capture_helpers)
target_link_libraries(testupdatetotalshift capture_helpers)
target_link_libraries(testcheckargs cli)
target_link_libraries(testgetbeaterror analysis )
target_link_libraries(testerf mymath )

target_compile_definitions(testcalculateall
    PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME testerf COMMAND testerf)
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
