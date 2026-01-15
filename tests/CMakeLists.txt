
# tests/CMakeLists.txt
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
