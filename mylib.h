#include "crosscorint.h"
#include "myarr.h"
#include "parseargs.h"
#include <stdio.h>

void linreg(const double* xarr,
            const double* yarr,
            unsigned int ArrayLength,
            double* a,
            double* b,
            double* s);

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

void writefile(FILE* fp, int* array, unsigned int ArrayLength);
void writefileDouble(FILE* fp, double* array, unsigned int ArrayLength);

void calculateTotalFromFile(unsigned int n,
                            FILE* rawfile,
                            unsigned int ArrayLength,
                            double threshold);

void calculateTotal(unsigned int n,
                    double* maxpos,
                    unsigned int ArrayLength,
                    double threshold);

double getBeatError(const struct myarr* totalTick, double rate, int verbose);

int checkUIntArg(int name, unsigned int* value, char* optarg);

int checkFileArg(int name, FILE** fp, const char* optarg, const char* mode);
void fillReference(FILE* fpDefPeak,
                   struct myarr* reference,
                   unsigned int teeth);

int shiftHalf(unsigned int value, unsigned int ArrayLength);
int modSigned(int value, unsigned int ArrayLength);
int getInt(char* ptr);
double getDouble(char* ptr);
int getDoublesFromStdin(size_t max_count, double* arr);
int getIntsFromStdin(size_t max_count, int* arr);
