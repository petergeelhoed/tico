#pragma once

#include "myarr.h"

#include <stdio.h>

void calculateTotalFromFile(size_t count,
                            FILE* rawfile,
                            size_t ArrayLength,
                            double threshold,
                            double rate);

double getBeatError(const struct myarr* totalTick, double rate, int verbose);
