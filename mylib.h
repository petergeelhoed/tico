#include <alsa/asoundlib.h>
#include <fftw3.h>

void syncwrite(int* input, int NN, char* file);
void syncappend(int* input, int NN, FILE* file);
void rescale(int* total, int NN);

void linregd(const float* xarr,
             const float* yarr,
             int NN,
             double* a,
             double* b,
             double* s);
void linreg(
    const int* xarr, const int* yarr, int NN, double* a, double* b, double* s);

int getmaxpos(int* array, int NN);

void printheader(double b, int NN, int l, float beatError);

void printspaces(int maxpos,
                 int val,
                 char* spaces,
                 int mod,
                 int columns,
                 double a,
                 int cvalue);

void fit10secs(double* a,
               double* b,
               double* s,
               int i,
               int* maxvals,
               int* maxes,
               int cvalue,
               int npeaks);
void writefile(FILE* fp, int* array, int NN);
void calculateTotal(int n, int* maxpos, int NN, double threshold);
int getBeatError(int* totalTick, int NN, int verbose);
void crosscorint(int NN, int* array, int* ref, int* cross);
void* threadWrite(void* arr);
void* threadAppend(void* arr);
void writearray(int* arr, int NN, const char* file);
