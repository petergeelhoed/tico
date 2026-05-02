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

void fillReference(FILE* fpDefPeak,
                   struct myarr* reference,
                   unsigned int teeth);

void shiftBufferData(unsigned int* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals);

void processLogging(CapConfig* cfg,
                    AppResources* res,
                    unsigned int totalTime,
                    unsigned int writeInterval);

void fitAndPrint(unsigned int tickIndex,
                 unsigned int globalTickIndex,
                 struct myarr* cumulativeTick,
                 AppResources* res,
                 CapConfig* cfg,
                 unsigned int arrayLength,
                 unsigned int mod,
                 unsigned int currentColumns);

void rotateDerivativeWindow(AppResources* res,
                            size_t arrayLength,
                            int cumulativeShift);

int findMaxPosition(AppResources* res,
                    struct myarr* cumulativeTick,
                    unsigned int globalTickIndex,
                    unsigned int tickIndex,
                    size_t arrayLength,
                    CapConfig* cfg);

int updateTotalShiftIfNeeded(int cumulativeShift,
                             int peakOffset,
                             unsigned int globalTickIndex,
                             unsigned int tickIndex,
                             AppResources* res,
                             CapConfig* cfg);
