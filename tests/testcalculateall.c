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
        calculateTotalFromFile(7, rawfile, 16000, 3.0, 48000);
        fclose(rawfile);
    }
    return 0;
}
