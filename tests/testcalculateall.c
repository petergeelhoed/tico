#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "analysis.h"
#include "parseargs.h"

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR "."
#endif

// Magic number constants
#define CALC_TOTAL_INDEX 7
#define CALC_TOTAL_LENGTH 16000
#define CALC_TOTAL_FACTOR 3.0
#define CALC_TOTAL_RATE 48000

int main(void)
{
    FILE* rawfile = NULL;
    char path[BUFFER_SIZE];

    if (snprintf(path, sizeof(path), "%s/rawNan", TEST_DATA_DIR) >=
        (int)sizeof(path))
    {
        (void)fprintf(stderr, "fixture path too long\n");
        return EXIT_FAILURE;
    }

    if (checkFileArg('r', &rawfile, path, "r") != 0)
    {
        return EXIT_FAILURE;
    }

    if (rawfile)
    {
        calculateTotalFromFile(CALC_TOTAL_INDEX,
                               rawfile,
                               CALC_TOTAL_LENGTH,
                               CALC_TOTAL_FACTOR,
                               CALC_TOTAL_RATE);
        (void)fclose(rawfile);
    }
    return 0;
}
