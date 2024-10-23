#include "crosscorint.h"

void linregd(const float* xarr,
             const float* yarr,
             unsigned int NN,
             double* a,
             double* b,
             double* s);
void linreg(const int* xarr,
            const int* yarr,
            unsigned int NN,
            double* a,
            double* b,
            double* s);

unsigned int getmaxpos(int* array, unsigned int NN);

void printheader(double b, unsigned int NN, unsigned int l, double beatError);

void printspaces(int maxpos,
                 int val,
                 char* spaces,
                 int mod,
                 unsigned int columns,
                 double a,
                 int cvalue);

void fit10secs(double* a,
               double* b,
               double* s,
               unsigned int i,
               int* maxvals,
               int* maxes,
               int cvalue,
               unsigned int npeaks);
void writefile(FILE* fp, int* array, unsigned int NN);
void calculateTotal(unsigned int n,
                    int* maxpos,
                    unsigned int NN,
                    double threshold);
unsigned int getBeatError(int* totalTick, unsigned int NN, int verbose);
int checkUIntArg(int name, unsigned int* value, char* optarg);
int checkFileArg(int name, FILE** fp, char* optarg, char* mode);
