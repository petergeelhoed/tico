#pragma once

#include "myarr.h"
#include <stddef.h>

/**
 * @defgroup mymath Mymath Library
 * Math and matrix utilities.
 * @{
 */

/** @ingroup mymath
 * @brief Transposes a matrix in-place.
 *
 * @param arr The matrix to transpose, stored as a 1D array in row-major order.
 * @param Nrows The number of rows in the original matrix.
 * @param Ncols The number of columns in the original matrix.
 */
void transpone(double* arr, unsigned int Nrows, unsigned int Ncols);

/** @ingroup mymath
 * @brief Inverts a square matrix in-place.
 *
 * @param arr The square matrix to invert, stored as a 1D array in row-major
 order. The result will overwrite the original matrix.
 * @param Nrows The number of rows (and columns) in the square matrix.
 * @param Ncols The number of columns (and rows) in the square matrix.
 */
void invert(double* arr, unsigned int Nrows, unsigned int Ncols);

/** @ingroup mymath
 * @brief Multiplies two matrices and returns the result as a new matrix.
 *
 * @param matrix0 The first matrix, stored as a 1D array in row-major order.
 * @param Nrows The number of rows in the first matrix.
 * @param Ncols The number of columns in the first matrix.
 * @param matrix1 The second matrix, stored as a 1D array in row-major order.
 * @param Mrows The number of rows in the second matrix.
 * @param Mcols The number of columns in the second matrix.
 * @return A pointer to the resulting matrix, stored as a 1D array in row-major
 * order. The caller is responsible for freeing this memory.
 */
double* mulmat(const double* matrix0,
               unsigned int Nrows,
               unsigned int Ncols,
               const double* matrix1,
               unsigned int Mrows,
               unsigned int Mcols);

/** @ingroup mymath
 * @brief Performs linear regression on a dataset with optional weighting.
 *
 * @param coeffs An array of size 2 where the resulting coefficients (intercept
 * and slope) will be stored.
 * @param xmat The design matrix, stored as a 1D array in row-major order.
 * @param Nrows The number of rows in the design matrix.
 * @param Ncols The number of columns in the design matrix.
 * @param vec The response vector, stored as a 1D array.
 * @param weight An optional array of weights for each data point, stored as a
 * 1D array. If NULL, unweighted regression is performed.
 */

void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int Nrows,
               unsigned int Ncols,
               double* vec,
               const double* weight);

/** @ingroup mymath
 * @brief Fits a line to N peaks in the data, with outlier rejection.
 *
 * @param intercept Pointer to store the intercept.
 * @param slope Pointer to store the slope.
 * @param curPos The current position in the data.
 * @param maxvals Array of weights.
 * @param maxes Array of y values.
 * @param subpos Array of sub-positions.
 * @param npeaks Number of peaks to fit.
 * @param SDthreshold Standard deviation threshold for outlier rejection.
 */
void fitNpeaks(double* intercept,
               double* slope,
               unsigned int curPos,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               unsigned int npeaks,
               double SDthreshold);

/** @ingroup mymath
 * @brief Performs fast weighted linear regression.
 *
 * @param coeffs Array to store the resulting coefficients (intercept and
 * slope).
 * @param xmat The array of x values.
 * @param Npoints The number of points.
 * @param vec The array of y values.
 * @param weightArr The array of weights.
 */
void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int Npoints,
                const double* vec,
                const double* weightArr);

/** @ingroup mymath
 * @brief Returns 1 if two doubles are nearly equal, 0 otherwise.
 *
 * @param number0 First number.
 * @param number1 Second number.
 * @return 1 if nearly equal, 0 otherwise.
 */
int nearlyEqual(double number0, double number1);

/** @ingroup mymath
 * @brief Returns the index of the maximum value in an integer array.
 *
 * @param array The array to search.
 * @param ArrayLength The number of elements in the array.
 * @return The index of the maximum value.
 */
size_t getmaxpos(const int* array, size_t ArrayLength);

/** @ingroup mymath
 * @brief Performs simple linear regression on two arrays.
 *
 * @param xarr The array of x values.
 * @param yarr The array of y values.
 * @param ArrayLength The number of elements in the arrays.
 * @param intercept Pointer to store the intercept.
 * @param slope Pointer to store the slope.
 * @param stdev Pointer to store the standard deviation of the fit.
 */
void linreg(const double* xarr,
            const double* yarr,
            size_t ArrayLength,
            double* intercept,
            double* slope,
            double* stdev);

/** @ingroup mymath
 * @brief Shifts a value by half the array length, wrapping around.
 *
 * @param value The value to shift.
 * @param ArrayLength The length of the array.
 * @return The shifted value, wrapped around the array length.
 */
int shiftHalf(size_t value, size_t ArrayLength);

/** @ingroup mymath
 * @brief Computes the modulus of a signed integer value with respect to the
 * length of an array, ensuring a non-negative result.
 *
 * @param value The signed integer value to mod.
 * @param ArrayLength The length of the array.
 * @return The modulus of the value with respect to the array length, wrapped to
 * be non-negative.
 */
size_t modSigned(int value, size_t ArrayLength);

/** @} */
