#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "../src/mylib.c"
#include "parseargs.h"

int main(void)
{
    FILE* rawfile;
    checkFileArg('r', &rawfile, "raw_nan", "r");

    if (rawfile)
    {
        calculateTotalFromFile(7, rawfile, 16000, 3.0, 48000);
        fclose(rawfile);
    }
    exit(0);
}
