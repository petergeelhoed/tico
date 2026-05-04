
#pragma once

#include "myarr.h"
#include <stdio.h>

/**
 * @defgroup analysis Analysis Library
 * Functions for analysis operations.
 * @{
 */

/** @ingroup analysis */
void calculateTotalFromFile(size_t count,
                            FILE* rawfile,
                            size_t ArrayLength,
                            double threshold,
                            double rate);

/** @ingroup analysis */
double getBeatError(const struct myarr* totalTick, double rate, int verbose);

/** @} */
