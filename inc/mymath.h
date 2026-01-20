#include "myarr.h"
void transpone(double*, unsigned int Nrows, unsigned int Ncols);
void invert(double*, unsigned int Nrows, unsigned int Ncols);

double* mulmat(const double* matrix0,
               unsigned int Nrows,
               unsigned int Ncols,
               const double* matrix1,
               unsigned int Mrows,
               unsigned int Mcols);

void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int Nrows,
               unsigned int Ncols,
               double* vec,
               const double* weight);

void fitNpeaks(double* intercept,
               double* slope,
               unsigned int cur_pos,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               unsigned int npeaks,
               double SDthreshold);

void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int Npoints,
                const double* vec,
                const double* weight);

int nearly_equal(double number0, double number1);

unsigned int getmaxpos(const int* array, unsigned int ArrayLength);

void linreg(const double* xarr,
            const double* yarr,
            unsigned int ArrayLength,
            double* intercept,
            double* slope,
            double* stdev);

int shiftHalf(unsigned int value, unsigned int ArrayLength);
int modSigned(int value, unsigned int ArrayLength);
