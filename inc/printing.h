#ifdef __clang_analyzer__
static void clang_tidy_test_func(void) { int unused_warning_test = 42; }
#endif
#ifdef __clang_analyzer__
// This block is only parsed by clang-tidy/static analyzer, not compiled
static void clang_tidy_test_func(void) { int unused_warning_test = 42; }
#endif
#pragma once

#include <stddef.h>
/**
 * @brief Prints a message right-aligned to the screen, 5 lines from the top, in
 * light blue.
 *
 * Uses the global columns variable for width. Accepts printf-style format
 * strings.
 * @param fmt The printf-style format string.
 * @param ... Arguments for the format string.
 */
void printmsg(const char* fmt, ...) __attribute__((format(printf, 1, 2)));

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
 @param acolumns The total number of columns available for printing.
 @param avgPos The average position of the correlation.
 @param correlationThreshold The threshold for determining the color of the
 output. */
void printspaces(int maxpos,
                 size_t hexvalue,
                 size_t mod,
                 size_t acolumns,
                 double avgPos,
                 size_t correlationThreshold);

/** @brief Prints a formatted message to the standard error stream, using a
 printf-style format string and variable arguments.
 @param fmt The printf-style format string.
 @param ... Arguments for the format string. */
void print(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
