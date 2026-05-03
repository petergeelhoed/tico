#pragma once
#include "config.h"
#include "myarr.h"
#include "mysound.h"

#include "fftw3.h"

#include <stdint.h>

/** @brief Struct to hold all the resources needed for the application. */
typedef struct
{
    struct myarr* subpos;
    struct myarr* maxpos;
    struct myarr* maxvals;
    struct myarr* derivative;
    struct myarr* tmpder;
    struct myarr* reference;
    struct myarr** teethArray;
    fftw_complex* filterFFT;
    int16_t* audioBuffer16;
} AppResources;

/** @brief Allocates and initializes the resources needed for the application.
 *
 * @param arrayLength The length of the arrays to be allocated.
 * @param ticktockBuffer The size of the audio buffer to be allocated.
 * @param cfg Pointer to the configuration structure containing necessary
 * parameters.
 * @return An AppResources struct containing pointers to all allocated
 * resources.
 */
AppResources allocateResources(size_t arrayLength,
                               size_t ticktockBuffer,
                               CapConfig* cfg);

/** @brief Frees all allocated resources and performs any necessary cleanup.
  @param res Pointer to the AppResources struct containing the resources to be
  freed.
  @param cfg Pointer to the CapConfig struct containing configuration parameters
  (if needed for cleanup).
  @param ctx Pointer to the CaptureCtx struct containing context information (if
  needed for cleanup).
 */
void cleanupResources(AppResources* res, CapConfig* cfg, CaptureCtx* ctx);
