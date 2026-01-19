#pragma once
#include <stdio.h>

#include "mydefs.h"

typedef struct
{
    double rate;
    unsigned int bph, evalue, zoom, time, everyline, cvalue, verbose, fitN,
        teeth;
    double SDthreshold;
    char device[MAX_DEVICE_LENGTH];
    FILE *fpposition, *fpmaxcor, *fptotal, *fpDefPeak, *fpInput;
} CapConfig;

int checkFileArg(int name,
                 FILE** filePtr,
                 const char* opt_arg,
                 const char* mode);

int checkUIntArg(int name, unsigned int* value, char* opt_arg);

void parse_arguments(int argc, char* argv[], CapConfig* cfg);
int getInt(char* ptr);
double getDouble(char* ptr);
int getDoublesFromStdin(size_t max_count, double* arr);
int getIntsFromStdin(size_t max_count, int* arr);
