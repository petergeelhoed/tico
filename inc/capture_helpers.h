
#pragma once
void print_card_device_mapping(void);

#include "config.h"
#include "myarr.h"
#include "resources.h"

#include <stdio.h>

/** Helper functions for capture.c, including printing, data shifting, logging,
 * and fitting. */

/* @brief Prints the header information for the capture output, including beat
error, fitted rate, and elapsed time.
@param fittedRate The fitted rate in seconds per day.
@param everyline A flag indicating whether to print the header on every line or
at the top of the screen
@param beatError The beat error in milliseconds.
@param seconds The elapsed time in seconds. */

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds);

/** @brief Prints spaces to visually represent the position of the maximum
 * correlation in the capture output, based on the provided parameters.
 @param maxpos The position of the maximum correlation.
 @param hexvalue The hexadecimal value representing the correlation strength.
 @param mod The modulus used for calculating the width of the spaces.
 @param columns The total number of columns available for printing.
 @param avgPos The average position of the correlation.
 @param correlationThreshold The threshold for determining the color of the
 output. */
void printspaces(int maxpos,
                 size_t hexvalue,
                 size_t mod,
                 size_t columns,
                 double avgPos,
                 size_t correlationThreshold);

/** @brief Prints the final results of the capture process, including total tick
 count, position data, and correlation data, based on the provided parameters.
 @param cfg The configuration settings for the capture process.
 @param res The resources used during the capture process.
 @param ArrayLength The length of the data arrays.
 @param totalTickTock The total tick count accumulated during the capture
 process. */
void printFinals(CapConfig* cfg,
                 AppResources* res,
                 size_t ArrayLength,
                 size_t totalTickTock);

/** @brief Fills the reference array with data from the provided file pointer,
 * based on the specified number of teeth.
 @param fpDefPeak The file pointer to read the reference data from.
 @param reference The myarr structure to store the reference data.
 @param teeth The number of teeth to determine how much data to read. */
void fillReference(FILE* fpDefPeak, struct myarr* reference, size_t teeth);

/** @brief Shifts the buffer data for subpos, maxpos, and maxvals arrays by
 * moving the data to the left and updating the ticktock value accordingly.
 @param ticktock A pointer to the current tick count that will be updated after
 shifting the data.
 @param subpos The myarr structure containing the subpos data to be shifted.
 @param maxpos The myarr structure containing the maxpos data to be shifted.
 @param maxvals The myarr structure containing the maxvals data to be shifted.
 */
void shiftBufferData(size_t* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals);

/** @brief Processes the logging of capture data, including printing the final
 results and handling any necessary logging based on the provided configuration
 and resources.
 @param cfg The configuration settings for the capture process, which may
 include file pointers for logging and other settings.
 @param res The resources used during the capture process, which may include
 data arrays and other relevant information.
 @param totalTime The total time elapsed during the capture process, which may
 be used for logging or final output.
 @param writeInterval The interval at which to write log data, which may be used
 to determine when to log certain information based on the total time. */
void processLogging(CapConfig* cfg,
                    AppResources* res,
                    size_t totalTime,
                    size_t writeInterval);

/** @brief Fits a model to the capture data and prints the results, including
 * calculating the intercept and slope, printing the header information, and
 * visualizing the position of the maximum correlation based on the provided
 * parameters.
 @param tickIndex The index of the current tick being processed.
 @param globalTickIndex The global index of the tick across the entire capture
 process.
 @param cumulativeTick The myarr structure containing cumulative tick data for
 beat error calculation.
 @param res The resources used during the capture process, which may include
 data arrays and other relevant information for fitting and printing.
 @param cfg The configuration settings for the capture process, which may
 include parameters for fitting and printing.
 @param arrayLength The length of the data arrays used in fitting and printing.
 @param mod The modulus used for calculating positions in printing.
 @param currentColumns The current number of columns available for printing,
 which may be used to determine how to visualize the output. */
void fitAndPrint(size_t tickIndex,
                 size_t globalTickIndex,
                 struct myarr* cumulativeTick,
                 AppResources* res,
                 CapConfig* cfg,
                 size_t arrayLength,
                 size_t mod,
                 size_t currentColumns);

/** @brief Rotates the derivative window by shifting the data in the derivative
 array based on the cumulative shift and array length, effectively updating the
 window for the next iteration of processing.
 @param res The resources used during the capture process, which may include
 data arrays and other relevant information for rotating the derivative window.
 @param arrayLength The length of the data arrays used in the derivative window,
 which is necessary for calculating the correct positions for rotation.
 @param cumulativeShift The total shift that has been accumulated, which will be
 used to determine how much to rotate the derivative window. */
void rotateDerivativeWindow(AppResources* res,
                            size_t arrayLength,
                            int cumulativeShift);

/** @brief Finds the maximum position of the correlation in the capture data by
 * performing an FFT fit and shifting the result based on the array length,
 using the provided resources and configuration settings.
 @param res The resources used during the capture process, which may include
 data arrays and other relevant information for finding the maximum position.
 @param cumulativeTick The myarr structure containing cumulative tick data for
 beat error calculation, which may be used in the process of finding the maximum
 position.
 @param globalTickIndex The global index of the tick across the entire capture
 process, which may be used to determine the current position in the data for
 finding the maximum position.
 @param tickIndex The index of the current tick being processed, which may be
 used to access specific data for finding the maximum position.
 @param arrayLength The length of the data arrays used in finding the maximum
 position, which is necessary for performing the FFT fit and shifting the result
 correctly.
 @param cfg The configuration settings for the capture process, which may
 include parameters for fitting and finding the maximum position, such as the
 number of peaks to fit and the standard deviation threshold. */
int findMaxPosition(AppResources* res,
                    struct myarr* cumulativeTick,
                    unsigned int globalTickIndex,
                    unsigned int tickIndex,
                    size_t arrayLength,
                    CapConfig* cfg);

/** @brief Updates the total shift if needed based on the cumulative shift, peak
 offset, and other parameters, using the provided resources and configuration
 settings.
 @param cumulativeShift The total shift that has been accumulated, which may be
 used to determine if an update to the total shift is needed based on the peak
 offset and other parameters.
 @param peakOffset The offset of the peak position, which may be used to
 determine if an update to the total shift is needed based on the cumulative
 shift and other parameters.
 @param globalTickIndex The global index of the tick across the entire capture
 process, which may be used to determine the current position in the data for
 updating the total shift if needed.
 @param tickIndex The index of the current tick being processed, which may be
 used to access specific data for updating the total shift if needed.
 @param res The resources used during the capture process, which may include
 data arrays and other relevant information for updating the total shift if
 needed.
 @param cfg The configuration settings for the capture process, which may
 include parameters for fitting and updating the total shift, such as thresholds
 and other relevant settings. */
int updateTotalShiftIfNeeded(int cumulativeShift,
                             int peakOffset,
                             size_t globalTickIndex,
                             size_t tickIndex,
                             AppResources* res,
                             CapConfig* cfg);
