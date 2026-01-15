#pragma once
#include <stdio.h>

#include "mydefs.h"

typedef struct
{
    double rate;
    unsigned int bph, evalue, zoom, time, everyline, cvalue, verbose, fitN,
        teeth;
    double SDthreshold;
    char* device;
    FILE *fpposition, *fpmaxcor, *fptotal, *fpDefPeak, *fpInput;
} CapConfig;

void parse_arguments(int argc, char* argv[], CapConfig* cfg);
