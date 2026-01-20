#pragma once
#include "config.h"

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
