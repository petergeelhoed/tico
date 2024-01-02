#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>


void syncwrite(int * input, int NN, char *file);
void syncappend(int * input, int NN, FILE* file);
snd_pcm_t * initAudio(snd_pcm_format_t format, char* device, unsigned int rate);
void rescale(int* total, int NN);

fftw_complex * makeFilter(int evalue, int NN);

int fftfit(int *mean, int *total, int *base, int *val, const fftw_complex *filterFFT, int NN,int verb);

void linreg(const int *xarr, const int *yarr, int NN, double *a, double *b, double *s);

void applyFilter50(int* input, int NN, double* out);
void applyFilter(int* input, int NN, fftw_complex* filterFFT, double* out);
fftw_complex* convolute(int NN, int *array, const fftw_complex *filter);
void remove50hz(int NN, int *array, int rate);

void normalise(int NN, fftw_complex *in);

int getmaxfftw(fftw_complex* array, int NN);
int getmaxpos(int * array, int NN);
void printspaces(int maxpos,int val, char* spaces,int mod,int columns, double a,double b,int NN, int cvalue, float beatError);

int readBuffer( snd_pcm_t *capture_handle, int NN, char *buffer, int *derivative);
void readBufferRaw( snd_pcm_t *capture_handle, int NN, char *buffer, int *in);

void fit10secs(double *a, double *b, double *s, int i,int* maxvals,int *maxes, int cvalue, int npeaks);
void writefiles(FILE* fptotal, FILE* rawfile, int* totaltick, int *maxpos, int n, int NN);
void calculateTotal(int n, int* maxpos,int NN, double threshold);
int readShiftedBuffer(int* derivative, snd_pcm_t *capture_handle, int NN, char* buffer, int maxpos, int shift, int* totalshift, int lowerBound, int upperBound);
fftw_complex* crosscor(int NN, fftw_complex* array, fftw_complex* ref);
int getBeatError(int* totalTick, int NN,int verbose);
void crosscorint(int NN, int* array, int* ref, int* cross);
void *threadWrite(void* arr);
void *threadAppend(void* arr);
void writearray(int* arr, int NN, const char* file);
void writearraydouble(double* arr, int NN, const char* file);
void writefftw(fftw_complex * arr, int NN, const char* file);
