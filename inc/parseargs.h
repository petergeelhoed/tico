#pragma once
#include "config.h"

/**
 * Parses command line arguments and fills the CapConfig struct.
 * Also provides utility functions for checking and parsing arguments.
 */

/** @brief Checks if the provided argument is a valid file and opens it.
 *
 * @param name The name of the argument (for error messages).
 * @param filePtr Pointer to a FILE* that will be set to the opened file.
 * @param optArg The argument string to check (should be a filename).
 * @param mode The mode to open the file (e.g., "r" for read, "w" for write).
 * @return 0 on success, non-zero on failure.
 */
int checkFileArg(int name,
                 FILE** filePtr,
                 const char* optArg,
                 const char* mode);

/** @brief Checks if the provided argument is a valid integer and parses it.
 *
 * @param name The name of the argument (for error messages).
 * @param value Pointer to an int that will be set to the parsed value.
 * @param optArg The argument string to check (should be an integer).
 * @return 0 on success, non-zero on failure.
 */
int checkUIntArg(int name, unsigned int* value, const char* optArg);

/** @brief Parses command line arguments and fills the CapConfig struct.
 *
 * @param argc The number of command line arguments.
 * @param argv The array of command line argument strings.
 * @param cfg Pointer to a CapConfig struct that will be filled with the parsed
 * configuration.
 */
void parseArguments(int argc, char* argv[], CapConfig* cfg);

/** @brief Parses an integer from a string.
 *
 * @param ptr The string to parse.
 * @return The parsed integer value.
 */
int getInt(char* ptr);

/** @brief Parses a double from a string.
 *
 * @param ptr The string to parse.
 * @return The parsed double value.
 */
double getDouble(char* ptr);

/** @brief Reads doubles from standard input until EOF or maxCount is reached.
 *
 * @param maxCount The maximum number of doubles to read.
 * @param arr Pointer to an array where the read doubles will be stored.
 * @return The number of doubles successfully read.
 */
int getDoublesFromStdin(size_t maxCount, double* arr);

/** @brief Reads integers from standard input until EOF or maxCount is reached.
 *
 * @param maxCount The maximum number of integers to read.
 * @param arr Pointer to an array where the read integers will be stored.
 * @return The number of integers successfully read.
 */
int getIntsFromStdin(size_t maxCount, int* arr);
