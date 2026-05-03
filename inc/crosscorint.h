#pragma once
#include <stddef.h>

/** @brief Computes the cross-correlation of two integer arrays.
 *
 * @param ArrayLength The length of the input arrays.
 * @param array The first input array.
 * @param ref The second input array (reference).
 * @param cross The output array where the cross-correlation results will be
 * stored.
 */
void crosscorint(size_t ArrayLength,
                 const int* array,
                 const int* ref,
                 int* cross);
