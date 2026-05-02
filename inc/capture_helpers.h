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
                 size_t ArrayLength,
                 size_t totalTickTock);

void fillReference(FILE* fpDefPeak, struct myarr* reference, size_t teeth);

void shiftBufferData(size_t* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals);

void processLogging(CapConfig* cfg,
                    AppResources* res,
                    size_t totalTime,
                    size_t writeInterval);

void fitAndPrint(size_t tickIndex,
                 size_t globalTickIndex,
                 struct myarr* cumulativeTick,
                 AppResources* res,
                 CapConfig* cfg,
                 size_t arrayLength,
                 size_t mod,
                 size_t currentColumns);

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
                             size_t globalTickIndex,
                             size_t tickIndex,
                             AppResources* res,
                             CapConfig* cfg);
