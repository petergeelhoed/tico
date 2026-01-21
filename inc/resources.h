#pragma once
#include "myarr.h"

#include "fftw3.h"

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
    char* audioBuffer;
} AppResources;
