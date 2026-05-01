
#include "parseargs.h"
#include <stdio.h>
#include <stdlib.h>

// Magic number constants
#define TEST_VALID_INPUT "42"
#define TEST_VALID_VALUE 42
#define TEST_ZERO_INPUT "0"
#define TEST_NEGATIVE_INPUT "-3"
#define TEST_NONNUMERIC_INPUT "abc"
#define TEST_BAD_FILENAME "-bad"
#define TEST_MISSING_FILENAME "definitely_missing_file_for_test"
#define RETURN_VALID_INPUT_FAIL 1
#define RETURN_ZERO_ACCEPTED 2
#define RETURN_NEGATIVE_ACCEPTED 3
#define RETURN_NONNUMERIC_ACCEPTED 4
#define RETURN_BAD_FILENAME_ACCEPTED 5
#define RETURN_MISSING_FILE_ACCEPTED 6

int main(void)
{
    unsigned int value = 0;

    if (checkUIntArg('x', &value, TEST_VALID_INPUT) != 0 ||
        value != TEST_VALID_VALUE)
    {
        (void)fprintf(stderr, "checkUIntArg valid input failed\n");
        return RETURN_VALID_INPUT_FAIL;
    }
    if (checkUIntArg('x', &value, TEST_ZERO_INPUT) == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted zero\n");
        return RETURN_ZERO_ACCEPTED;
    }
    if (checkUIntArg('x', &value, TEST_NEGATIVE_INPUT) == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted negative value\n");
        return RETURN_NEGATIVE_ACCEPTED;
    }
    if (checkUIntArg('x', &value, TEST_NONNUMERIC_INPUT) == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted non-numeric value\n");
        return RETURN_NONNUMERIC_ACCEPTED;
    }

    FILE* filePtr = NULL;
    if (checkFileArg('r', &filePtr, TEST_BAD_FILENAME, "r") == 0)
    {
        (void)fprintf(stderr, "checkFileArg accepted option-like filename\n");
        return RETURN_BAD_FILENAME_ACCEPTED;
    }

    if (checkFileArg('r', &filePtr, TEST_MISSING_FILENAME, "r") == 0)
    {
        (void)fprintf(stderr, "checkFileArg accepted missing file\n");
        return RETURN_MISSING_FILE_ACCEPTED;
    }

    return 0;
}
