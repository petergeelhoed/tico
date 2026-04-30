#pragma once
#include "config.h"

int checkFileArg(int name,
                 FILE** filePtr,
                 const char* optArg,
                 const char* mode);

int checkUIntArg(int name, unsigned int* value, char* optArg);

void parseArguments(int argc, char* argv[], CapConfig* cfg);
int getInt(char* ptr);
double getDouble(char* ptr);
int getDoublesFromStdin(size_t maxCount, double* arr);
int getIntsFromStdin(size_t maxCount, int* arr);
