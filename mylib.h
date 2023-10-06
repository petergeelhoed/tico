#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

snd_pcm_t * initAudio(snd_pcm_format_t format, char* device, unsigned int rate);

fftw_complex * makeFilter(int evalue, int buffer_frames);

int fftfit(int *mean, int *total, FILE* rawfile, int *base, int *val, const fftw_complex *filterFFT, int buffer_frames);


