
#pragma once
#include <stdlib.h>

/**
 * @defgroup myarr Myarr Library
 * Functions for array management.
 * @{
 */

/** @brief A struct to hold an array of integers or an array of doubles, along
 * with the length of the arrays. */
struct myarr
{
    int* arr;
    double* arrd;
    size_t ArrayLength;
};

/** @ingroup myarr @brief Frees the memory allocated for the arrays in a myarr
 * struct.
 * @param arr A pointer to the myarr struct whose arrays are to be freed. */
void freemyarr(struct myarr* arr);

/** @ingroup myarr @brief Creates a myarr struct with an array of integers of
 * the specified length.
 * @param ArrayLength The length of the integer array to be created.
 * @return A pointer to the newly created myarr struct with the integer array.
 */
struct myarr* makemyarr(size_t ArrayLength);

/** @ingroup myarr @brief Creates a myarr struct with an array of doubles of the
 * specified length.
 * @param ArrayLength The length of the double array to be created.
 * @return A pointer to the newly created myarr struct with the double array. */
struct myarr* makemyarrd(size_t ArrayLength);

/** @} */
