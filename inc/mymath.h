#include "myarr.h"

/** @brief Transposes a matrix in-place.
 *
 * @param arr The matrix to transpose, stored as a 1D array in row-major order.
 * @param Nrows The number of rows in the original matrix.
 * @param Ncols The number of columns in the original matrix.
 */
void transpone(double* arr, unsigned int Nrows, unsigned int Ncols);

/** @brief Inverts a square matrix in-place.

 * @param arr The square matrix to invert, stored as a 1D array in row-major
 order. The result will overwrite the original matrix.
 * @param Nrows The number of rows (and columns) in the square matrix.
 * @param Ncols The number of columns (and rows) in the square matrix.
 */
void invert(double* arr, unsigned int Nrows, unsigned int Ncols);

/** @brief Multiplies two matrices and returns the result as a new matrix.
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

/** @brief Performs linear regression on a dataset with optional weighting.
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

/** @brief Fits a linear model to a specified number of peaks in the data.
 *
 * @param intercept A pointer where the resulting intercept of the fitted line
 * will be stored.
 * @param slope A pointer where the resulting slope of the fitted line will be
 * stored.
 * @param curPos The current position in the data from which to start fitting.
 * @param maxvals An array containing the maximum values of the peaks.
 * @param maxes An array containing the positions of the peaks.
 * @param subpos An array containing the positions of the sub-peaks or other
 * relevant positions.
 * @param npeaks The number of peaks to fit.
 * @param SDthreshold A threshold for standard deviation to determine which
 * peaks to include in the fit.
 */
void fitNpeaks(double* intercept,
               double* slope,
               unsigned int curPos,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               unsigned int npeaks,
               double SDthreshold);

/** @brief A faster implementation of linear regression that may be less
 * accurate than matlinreg.
 *
 * @param coeffs An array of size 2 where the resulting coefficients (intercept
 * and slope) will be stored.
 * @param xmat The design matrix, stored as a 1D array in row-major order.
 * @param Npoints The number of data points (rows) in the design matrix.
 * @param vec The response vector, stored as a 1D array.
 * @param weight An optional array of weights for each data point, stored as a
 * 1D array. If NULL, unweighted regression is performed.
 */
void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int Npoints,
                const double* vec,
                const double* weight);

/** @brief Checks if two double precision floating-point numbers are nearly
 * equal, accounting for potential precision issues.
 *
 * @param number0 The first number to compare.
 * @param number1 The second number to compare.
 * @return 1 if the numbers are nearly equal, 0 otherwise.
 */
int nearlyEqual(double number0, double number1);

/** @brief Finds the position of the maximum value in an array of integers.
 *
 * @param array The array to search through.
 * @param ArrayLength The length of the array.
 * @return The index of the maximum value in the array.
 */
size_t getmaxpos(const int* array, size_t ArrayLength);

/** @brief Performs linear regression on two arrays of doubles and computes the
 * intercept, slope, and standard deviation of the fit.
 *
 * @param xarr The array of x-values.
 * @param yarr The array of y-values.
 * @param ArrayLength The length of the x and y arrays.
 * @param intercept A pointer where the resulting intercept of the fitted line
 * will be stored.
 * @param slope A pointer where the resulting slope of the fitted line will be
 * stored.
 * @param stdev A pointer where the resulting standard deviation of the fit will
 * be stored.
 */
void linreg(const double* xarr,
            const double* yarr,
            size_t ArrayLength,
            double* intercept,
            double* slope,
            double* stdev);

/** @brief Shifts a value by half the length of an array, wrapping around.
 *
 * @param value The value to shift.
 * @param ArrayLength The length of the array.
 * @return The shifted value, wrapped around the array length.
 */
int shiftHalf(size_t value, size_t ArrayLength);

/** @brief Computes the modulus of a signed integer value with respect to the
 * length of an array, ensuring a non-negative result.
 *
 * @param value The signed integer value to mod.
 * @param ArrayLength The length of the array.
 * @return The modulus of the value with respect to the array length, wrapped to
 * be non-negative.
 */
size_t modSigned(int value, size_t ArrayLength);
