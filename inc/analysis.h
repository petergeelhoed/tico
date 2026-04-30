#pragma once

#include "myarr.h"

#include <stdio.h>

void calculateTotalFromFile(unsigned int count,
                            FILE* rawfile,
                            unsigned int ArrayLength,
                            double threshold,
                            double rate);

double getBeatError(const struct myarr* totalTick, double rate, int verbose);
