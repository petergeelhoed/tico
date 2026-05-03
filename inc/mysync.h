#include "myarr.h"

#include <stdio.h>

/** @brief Synchronization functions for myarr struct and file operations. */

/** @brief Waits for the file to be closed before proceeding.
  @param file Pointer to the file to wait on. */
void waitClose(FILE* file);

/** @brief Waits for all pending operations to complete before proceeding.
  checks locks for ctr_mutex and waits on ctr_zero until count is zero, then
  releases lock. */
void wait(void);

/** @brief Writes the input array to the specified file in a separate thread.
  @param input Pointer to the input array to write.
  @param ArrayLength Length of the input array.
  @param file Name of the file to write to.
  @return 0 on success, -1 on failure. */
int syncwrite(int* input, size_t ArrayLength, const char* file);

/** @brief Prints the current time of day to the specified file with
  synchronization.
  @param out Pointer to the file to write to. */
void printTOD(FILE* out);

/** @brief Appends the input myarr struct to the specified file in a separate
  thread.
  @param input Pointer to the input myarr struct to append.
  @param file Pointer to the file to append to. */
void syncAppendMyarr(struct myarr* input, FILE* file);

/** @brief Worker function for appending a myarr struct to a file in a separate
  thread.
  @param inStruct Pointer to the append_task struct containing the myarr and
  file information.
  @return NULL on completion. */
void* threadAppendMyarr(void* inStruct);
