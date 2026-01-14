#include "myarr.h"
#include <fftw3.h>

fftw_complex* makeFilter(unsigned int evalue, unsigned int ArrayLength);

unsigned int fftfit(struct myarr input,
                    int* total,
                    const int* base,
                    double* corvalue,
                    fftw_complex* filterFFT,
                    int verb,
                    double* subpos);
fftw_complex* convolute(struct myarr array, fftw_complex* filter);

unsigned int getmaxfftw(fftw_complex* array, unsigned int ArrayLength);
fftw_complex* crosscor(unsigned int ArrayLength, fftw_complex* array, fftw_complex* ref);
void writefftw(fftw_complex* arr, unsigned int ArrayLength, const char* file);

void remove50hz(unsigned int ArrayLength, int* array, unsigned int rate);
void normalise(unsigned int ArrayLength, fftw_complex* in);
void rescale(int* total, unsigned int ArrayLength);
int getshift(struct myarr xarr, struct myarr yarr);
