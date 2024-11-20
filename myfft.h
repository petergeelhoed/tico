#include <fftw3.h>
#include "myarr.h"

fftw_complex* makeFilter(unsigned int evalue, unsigned int NN);

unsigned fftfit(const struct myarr input,
                int* total,
                int* base,
                double* val,
                fftw_complex* filterFFT,
                unsigned int NN,
                int verb);
void applyFilter(int* input,
                 unsigned int NN,
                 fftw_complex* filterFFT,
                 double* out);
fftw_complex* convolute(unsigned int NN, int* array, fftw_complex* filter);

unsigned int getmaxfftw(fftw_complex* array, unsigned int NN);
fftw_complex* crosscor(unsigned int NN, fftw_complex* array, fftw_complex* ref);
void writefftw(fftw_complex* arr, unsigned int NN, const char* file);

void remove50hz(unsigned int NN, int* array, unsigned int rate);
void normalise(unsigned int NN, fftw_complex* in);
void rescale(int* total, unsigned int NN);
