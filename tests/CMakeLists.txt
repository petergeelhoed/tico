
# tests/CMakeLists.txt

set(NO_RUN_TESTS
    testmicamplification
    testfft
)
set(TEST_ERROUT_TARGETS 
    testcalculateall
   )

set(TEST_EXPOUT_TARGETS 
    testcrosscorint
    testinvert
   )
set(TEST_TARGETS 
    testsyncwrite
    testssecshift
    testsyncappend
    testgetdoubles
    testmysync
    testfilter
    testremove50
    testlinreg 
    testmulmat
    testmatlinreg
    testmatlinregquad 
    testtranspone
    testshiftbufferdata
    testrotatederivativewindow
    testupdatetotalshift
    testcheckargs
    testgetbeaterror
    testgetinput
    testerf
    testcomparedoubles
)
foreach(TEST
        ${TEST_TARGETS}
        ${NO_RUN_TESTS} 
        ${TEST_EXPOUT_TARGETS}
        ${TEST_ERROUT_TARGETS}
        )
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
target_link_libraries(testcomparedoubles mymath )
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
target_link_libraries(testssecshift mymath )

target_compile_definitions(testcalculateall
    PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

set( TESTS
        ${TEST_TARGETS}
        ${TEST_EXPOUT_TARGETS} 
        ${TEST_ERROUT_TARGETS} 
        )

foreach(TST ${TEST_TARGETS})
    add_test(NAME ${TST} COMMAND $<TARGET_FILE:${TST}>)
endforeach()

foreach(TST ${TEST_ERROUT_TARGETS})
    add_test(
      NAME ${TST}
      COMMAND /bin/sh -c
              "set -e; \
               \"$<TARGET_FILE:${TST}>\" 2> \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\"; \
               \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_${TST}.txt\""
    )
endforeach()

foreach(TST ${TEST_EXPOUT_TARGETS})
    add_test(
      NAME ${TST}
      COMMAND /bin/sh -c
              "set -e; \
               \"$<TARGET_FILE:${TST}>\" > \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\"; \
               \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_${TST}.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_${TST}.txt\""
    )
endforeach()
