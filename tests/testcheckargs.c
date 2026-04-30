#include <stdio.h>
#include <stdlib.h>

#include "parseargs.h"

int main(void)
{
    unsigned int value = 0;

    if (checkUIntArg('x', &value, "42") != 0 || value != 42)
    {
        (void)fprintf(stderr, "checkUIntArg valid input failed\n");
        return 1;
    }
    if (checkUIntArg('x', &value, "0") == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted zero\n");
        return 2;
    }
    if (checkUIntArg('x', &value, "-3") == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted negative value\n");
        return 3;
    }
    if (checkUIntArg('x', &value, "abc") == 0)
    {
        (void)fprintf(stderr, "checkUIntArg accepted non-numeric value\n");
        return 4;
    }

    FILE* filePtr = NULL;
    if (checkFileArg('r', &filePtr, "-bad", "r") == 0)
    {
        (void)fprintf(stderr, "checkFileArg accepted option-like filename\n");
        return 5;
    }

    if (checkFileArg('r', &filePtr, "definitely_missing_file_for_test", "r") ==
        0)
    {
        (void)fprintf(stderr, "checkFileArg accepted missing file\n");
        if (filePtr != NULL)
        {
            (void)fclose(filePtr);
        }
        return 6;
    }

    FILE* seed = fopen("tmp_checkfilearg.txt", "w");
    if (seed == NULL)
    {
        (void)fprintf(stderr, "cannot create temp file for checkFileArg test\n");
        return 7;
    }
    (void)fputs("ok\n", seed);
    (void)fclose(seed);

    filePtr = NULL;
    if (checkFileArg('r', &filePtr, "tmp_checkfilearg.txt", "r") != 0 ||
        filePtr == NULL)
    {
        (void)fprintf(stderr, "checkFileArg failed on valid path\n");
        (void)remove("tmp_checkfilearg.txt");
        return 8;
    }

    (void)fclose(filePtr);
    (void)remove("tmp_checkfilearg.txt");
    return 0;
}
