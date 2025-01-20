#include "crosscorint.h"
#include "myarr.h"

void linregd(const double* xarr,
             const double* yarr,
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

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double tps);

void printspaces(int maxpos,
                 double hexvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double a,
                 unsigned int cvalue);

void fit10secs(double* a,
               double* b,
               double* s,
               const unsigned int i,
               const int* maxvals,
               const int* maxes,
               const int cvalue,
               const unsigned int npeaks);

void writefile(FILE* fp, int* array, unsigned int NN);
void writefileDouble(FILE* fp, double* array, unsigned int NN);

void calculateTotalFromFile(unsigned int n,
                            FILE* rawfile,
                            unsigned int NN,
                            double threshold);

void calculateTotal(unsigned int n,
                    int* maxpos,
                    unsigned int NN,
                    double threshold);

double
getBeatError(const struct myarr* totalTick, unsigned int rate, int verbose);

int checkUIntArg(int name, unsigned int* value, char* optarg);

int checkFileArg(int name, FILE** fp, char* optarg, char* mode);
void fillReference(FILE* fpDefPeak,
                   struct myarr* reference,
                   unsigned int teeth);

void checkAndFlip(struct myarr* totaltick,
                  struct myarr* reference,
                  unsigned int verbose);
int shiftHalf(unsigned int value, unsigned int NN);
int modSigned(int value, unsigned int NN);
