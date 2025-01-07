#include "myarr.h"
void transpone(double*, unsigned int N, unsigned int M);
void invert(double*, unsigned int N, unsigned int M);

double* mulmat(double*,
               unsigned int N,
               unsigned int M,
               double*,
               unsigned int S,
               unsigned int T);
void matlinreg(double coeffs[2],
               double* xmat,
               unsigned int N,
               unsigned int M,
               double* y,
               double* w);

void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const unsigned int npeaks);
