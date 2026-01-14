#include "myarr.h"
void transpone(double*, unsigned int Nrows, unsigned int Ncols);
void invert(double*, unsigned int Nrows, unsigned int Ncols);

double* mulmat(const double* matrix,
               unsigned int Nrows,
               unsigned int Ncols,
               const double* vector,
               unsigned int S,
               unsigned int T);
void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int Nrows,
               unsigned int Ncols,
               double* vec,
               const double* weight);

void fitNpeaks(double* par_a,
               double* par_b,
               unsigned int i,
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
