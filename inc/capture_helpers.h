#pragma once

#include "config.h"
#include "myarr.h"
#include "resources.h"

#include <stdio.h>

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds);

void printspaces(int maxpos,
                 double hexvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double avgPos,
                 unsigned int correlationThreshold);

void printFinals(CapConfig* cfg,
                 AppResources* res,
                 unsigned int ArrayLength,
                 unsigned int totalTickTock);

void fillReference(FILE* fpDefPeak, struct myarr* reference, unsigned int teeth);

void shiftBufferData(unsigned int* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals);
