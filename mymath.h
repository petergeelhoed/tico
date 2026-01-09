#include "myarr.h"
void transpone(double*, unsigned int N, unsigned int M);
void invert(double*, unsigned int N, unsigned int M);

double* mulmat(const double* matrix,
               unsigned int N,
               unsigned int M,
               const double* vector,
               unsigned int S,
               unsigned int T);
void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int N,
               unsigned int M,
               double* vec,
               const double* weight);

void fitNpeaks(double* a,
               double* b,
               unsigned int i,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               unsigned int npeaks);

void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int N,
                const double* vec,
                const double* weight);
