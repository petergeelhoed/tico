
#pragma once

#include "myarr.h"
#include <stdio.h>

/**
 * @defgroup analysis Analysis Library
 * Functions for analysis operations.
 * @{
 */

/**
 * @ingroup analysis
 * @brief Computes total values from a raw input file.
 * @param count Number of records to process.
 * @param rawfile Input file containing raw data.
 * @param ArrayLength Length of each data array.
 * @param threshold Threshold used to select values.
 * @param rate Sample rate associated with the input data.
 */
void calculateTotalFromFile(size_t count,
                            FILE* rawfile,
                            size_t ArrayLength,
                            double threshold,
                            double rate);

/**
 * @ingroup analysis
 * @brief Calculates beat error from the total tick data.
 * @param totalTick Array containing total tick values.
 * @param rate Sample rate used to interpret the tick data.
 * @param verbose Nonzero to print diagnostic details.
 * @return The calculated beat error.
 */
double getBeatError(const struct myarr* totalTick, double rate, int verbose);

/** @} */
