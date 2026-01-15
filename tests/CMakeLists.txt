
# tests/CMakeLists.txt

add_test(
  NAME testinvert
  COMMAND /bin/sh -c
          "set -e; \
           \"$<TARGET_FILE:testinvert>\" > \"${CMAKE_CURRENT_BINARY_DIR}/actual_testinvert.txt\"; \
           \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_testinvert.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_testinvert.txt\""
)

