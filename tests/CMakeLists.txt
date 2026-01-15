
# tests/CMakeLists.txt

add_test(
  NAME testinvert
  COMMAND /bin/sh -c
          "set -e; \
           \"$<TARGET_FILE:noon>\" < \"${CMAKE_CURRENT_BINARY_DIR}/input_testnoon.txt\" > \"${CMAKE_CURRENT_BINARY_DIR}/actual_testnoon.txt\"; \
           \"${CMAKE_COMMAND}\" -E compare_files \"${CMAKE_CURRENT_BINARY_DIR}/actual_testnoon.txt\" \"${CMAKE_CURRENT_SOURCE_DIR}/expected_testnoon.txt\""
)

