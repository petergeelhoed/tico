#include "config.h"
#include "crosscorint.h"
#include "myarr.h"
#include <stdio.h>

void linreg(const double* xarr,
            const double* yarr,
            unsigned int ArrayLength,
            double* intercept,
            double* slope,
            double* stdev);

unsigned int getmaxpos(const int* array, unsigned int ArrayLength);

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds);

void printspaces(int maxpos,
                 double hexvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double avg_pos,
                 unsigned int correlationThreshold);

void writefile(FILE* filePtr, int* array, unsigned int ArrayLength);
void writefileDouble(FILE* filePtr, double* array, unsigned int ArrayLength);

void calculateTotalFromFile(unsigned int count,
                            FILE* rawfile,
                            unsigned int ArrayLength,
                            double threshold,
                            double rate);

double getBeatError(const struct myarr* totalTick, double rate, int verbose);

void fillReference(FILE* fpDefPeak,
                   struct myarr* reference,
                   unsigned int teeth);

int shiftHalf(unsigned int value, unsigned int ArrayLength);
int modSigned(int value, unsigned int ArrayLength);
