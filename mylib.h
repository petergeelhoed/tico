#include "crosscorint.h"

void linregd(const float* xarr,
             const float* yarr,
             int NN,
             double* a,
             double* b,
             double* s);
void linreg(
    const int* xarr, const int* yarr, int NN, double* a, double* b, double* s);

unsigned int getmaxpos(int* array, unsigned int NN);

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
unsigned int getBeatError(int* totalTick, unsigned int NN, int verbose);
int checkUIntArg(int name, unsigned int* value, char* optarg);
int checkFileArg(int name, FILE*fp,  char* optarg, char* mode);
