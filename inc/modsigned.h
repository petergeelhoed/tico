#pragma once

#include <stddef.h>
/** * @brief Computes the modulus of a signed integer with respect to an array
 * length.
 *
 * This function takes a signed integer value and an array length, and returns
 * the modulus of the value with respect to the array length. The result is
 * always a non-negative integer that falls within the range of valid indices
 * for the array.
 *
 * @param value The signed integer value to be modded.
 * @param ArrayLength The length of the array, which must be greater than 0.
 * @return The modulus of the value with respect to the array length, as a
 * size_t.
 */
size_t modSigned(int value, size_t ArrayLength);
