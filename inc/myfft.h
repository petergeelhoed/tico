#pragma once

#include "myarr.h"
#include <fftw3.h>

/**
 * @defgroup myfft Myfft Library
 * FFT and filter utilities.
 * @{
 */

/** @ingroup myfft @brief Make a filter for the given evalue and array length.
 *
 * @param evalue The evalue to use for the filter.
 * @param ArrayLength The length of the array to filter.
 * @return A pointer to the filter in the frequency domain.
 */
fftw_complex* makeFilter(size_t evalue, size_t ArrayLength);

/** @ingroup myfft @brief Fit the input array to the filter and return the total
 * shift.
 *
 * @param input The input array to fit.
 * @param total A pointer to an integer where the total shift will be stored.
 * @param base A pointer to an integer where the base shift will be stored.
 * @param corvalue A pointer to a double where the correlation value will be
 * stored.
 * @param filterFFT The filter in the frequency domain.
 * @param verb An integer indicating whether to print verbose output (1 for yes,
 * 0 for no).
 * @param subpos A pointer to a double where the subpixel position will be
 * stored.
 * @return The total shift as a size_t.
 */
size_t fftfit(struct myarr input,
              int* total,
              const int* base,
              double* corvalue,
              fftw_complex* filterFFT,
              int verb,
              double* subpos);

/** @ingroup myfft @brief Convolute the input array with the given filter in the
 * frequency domain.
 *
 * @param array The input array to convolute.
 * @param filter The filter in the frequency domain to convolute with.
 * @return A pointer to the result of the convolution in the frequency domain.
 */
fftw_complex* convolute(struct myarr array, fftw_complex* filter);

/** @ingroup myfft @brief Get the maximum value from the given array in the
 * frequency domain.
 *
 * @param array The input array in the frequency domain to search for the
 * maximum value.
 * @param ArrayLength The length of the input array.
 * @return The maximum value found in the input array as a size_t.
 */
size_t getmaxfftw(fftw_complex* array, size_t ArrayLength);

/** @ingroup myfft
 * @brief Perform cross-correlation between the input array and the reference
 * array in the frequency domain.
 *
 * @param ArrayLength The length of the input and reference arrays.
 * @param array The input array in the frequency domain to cross-correlate.
 * @param ref The reference array in the frequency domain to cross-correlate
 * with.
 * @return A pointer to the result of the cross-correlation in the frequency
 * domain.
 */
fftw_complex* crosscor(size_t ArrayLength,
                       fftw_complex* array,
                       fftw_complex* ref);

/** @ingroup myfft
 * @brief Write the given array in the frequency domain to a file.
 *
 * @param arr The input array in the frequency domain to write to a file.
 * @param ArrayLength The length of the input array.
 * @param file The name of the file to write the array to.
 */
void writefftw(fftw_complex* arr, size_t ArrayLength, const char* file);

/** @ingroup myfft
 * @brief Remove the 50Hz component from the input array.
 *
 * @param ArrayLength The length of the input array.
 * @param array The input array to remove the 50Hz component from.
 * @param rate The sampling rate of the input array.
 */
void remove50hz(size_t ArrayLength, int* array, unsigned int rate);

/** @ingroup myfft
 * @brief Normalise the input array in the frequency domain.
 *
 * @param ArrayLength The length of the input array.
 * @param inData The input array in the frequency domain to normalise.
 */
void normalise(size_t ArrayLength, fftw_complex* inData);

/** @ingroup myfft
 * @brief Rescale the input array to fit within the range of a 16-bit signed
 * integer.
 *
 * @param total A pointer to an integer where the total shift will be stored.
 * @param ArrayLength The length of the input array.
 */
void rescale(int* total, size_t ArrayLength);

/** @ingroup myfft
 * @brief Get the shift between the two input arrays.
 *
 * @param xarr The first input array to compare.
 * @param yarr The second input array to compare.
 * @return The shift between the two input arrays as an integer.
 */
int getshift(struct myarr xarr, struct myarr yarr);

/** @} */
