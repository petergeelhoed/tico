#pragma once
#include "mydefs.h"
#include <alsa/asoundlib.h>
#include <stdio.h>

/** @brief Configuration parameters for the capture process. */
typedef struct
{
    /** Requested capture sample rate. */
    double rate;
    /** Target beats per hour. */
    unsigned int bph;
    /** Gaussian convolution parameter. */
    unsigned int evalue;
    /** Display zoom factor. */
    unsigned int zoom;
    /** Maximum capture duration in seconds. */
    unsigned int time;
    /** Number of ticks between per-line output updates. */
    unsigned int everyline;
    /** Cross-correlation limit. */
    double cvalue;
    /** Peak selected for verbose per-peak output. */
    unsigned int verbose;
    /** Number of samples used for fitting. */
    unsigned int fitN;
    /** Number of teeth represented by the reference data. */
    unsigned int teeth;
    /** Number of ticks accumulated between writes. */
    unsigned int writeInterval;
    /** Standard-deviation threshold for filtering. */
    double SDthreshold;
    /** ALSA capture-device name. */
    char device[MAX_DEVICE_LENGTH];
    /** Output file for position data. */
    FILE* fpposition;
    /** Output file for maximum-correlation data. */
    FILE* fpmaxcor;
    /** Output file for total data. */
    FILE* fptotal;
    /** Input file containing the default peak reference. */
    FILE* fpDefPeak;
    /** Optional input audio file. */
    FILE* fpInput;
    /** Active ALSA capture handle. */
    snd_pcm_t* captureHandle;
    /** Pointer to application-wide runtime state. */
    struct AppState* appState;

} CapConfig;
