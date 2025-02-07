#include "myarr.h"
#include <fftw3.h>

fftw_complex* makeFilter(unsigned int evalue, unsigned int NN);

unsigned fftfit(const struct myarr input,
                int* total,
                int* base,
                double* corvalue,
                fftw_complex* filterFFT,
                int verb,
                double* subpos);
void applyFilter(const struct myarr input,
                 fftw_complex* filterFFT,
                 double* out);
fftw_complex* convolute(const struct myarr array, fftw_complex* filter);

unsigned int getmaxfftw(fftw_complex* array, unsigned int NN);
fftw_complex* crosscor(unsigned int NN, fftw_complex* array, fftw_complex* ref);
void writefftw(fftw_complex* arr, unsigned int NN, const char* file);

void remove50hz(unsigned int NN, int* array, unsigned int rate);
void normalise(unsigned int NN, fftw_complex* in);
void rescale(int* total, unsigned int NN);
int getshift(const struct myarr x, const struct myarr y);
