#include "myarr.h"
#include <fftw3.h>

fftw_complex* makeFilter(size_t evalue, size_t ArrayLength);

size_t fftfit(struct myarr input,
              int* total,
              const int* base,
              double* corvalue,
              fftw_complex* filterFFT,
              int verb,
              double* subpos);
fftw_complex* convolute(struct myarr array, fftw_complex* filter);

size_t getmaxfftw(fftw_complex* array, size_t ArrayLength);
fftw_complex* crosscor(size_t ArrayLength,
                       fftw_complex* array,
                       fftw_complex* ref);
void writefftw(fftw_complex* arr, size_t ArrayLength, const char* file);

void remove50hz(size_t ArrayLength, int* array, unsigned int rate);
void normalise(size_t ArrayLength, fftw_complex* inData);
void rescale(int* total, size_t ArrayLength);
int getshift(struct myarr xarr, struct myarr yarr);
