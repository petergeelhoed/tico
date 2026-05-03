#pragma once
#include "mydefs.h"
#include <alsa/asoundlib.h>
#include <stdio.h>

/** brief struct to hold the configuration parameters for the capture process */
typedef struct
{
    double rate;
    unsigned int bph;
    unsigned int evalue;
    unsigned int zoom;
    unsigned int time;
    unsigned int everyline;
    unsigned int cvalue;
    unsigned int verbose;
    unsigned int fitN;
    unsigned int teeth;
    double SDthreshold;
    char device[MAX_DEVICE_LENGTH];
    FILE* fpposition;
    FILE* fpmaxcor;
    FILE* fptotal;
    FILE* fpDefPeak;
    FILE* fpInput;
    snd_pcm_t* captureHandle;

} CapConfig;
