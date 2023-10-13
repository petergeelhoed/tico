#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

snd_pcm_t * initAudio(snd_pcm_format_t format, char* device, unsigned int rate);

fftw_complex * makeFilter(int evalue, int NN);

int fftfit(int *mean, int *total, int *base, int *val, const fftw_complex *filterFFT, int NN);

void linreg(const int *xarr, const int *yarr, int NN, double *a, double *b, double *s);

void convolute(int NN, int *array, const fftw_complex *filter);
