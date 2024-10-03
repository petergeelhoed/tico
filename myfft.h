#include <alsa/asoundlib.h>
#include <fftw3.h>


fftw_complex* makeFilter(int evalue, int NN);

int fftfit(int* mean,
           int* total,
           int* base,
           int* val,
           fftw_complex* filterFFT,
           int NN,
           int verb);
void applyFilter(int* input, int NN, fftw_complex* filterFFT, double* out);
fftw_complex* convolute(int NN, int* array, fftw_complex* filter);

int getmaxfftw(fftw_complex* array, int NN);
fftw_complex* crosscor(int NN, fftw_complex* array, fftw_complex* ref);
void writefftw(fftw_complex* arr, int NN, const char* file);

void remove50hz(int NN, int* array, int rate);
void normalise(int NN, fftw_complex* in);
void rescale(int* total, int NN);
