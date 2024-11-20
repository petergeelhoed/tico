void transpone(double*, unsigned int N, unsigned int M);
void invert(double*, unsigned int N, unsigned int M);

double* mulmat(double*,
               unsigned int N,
               unsigned int M,
               double*,
               unsigned int S,
               unsigned int T);
double* matlinreg(
               double* xmat,
               unsigned int N,
               unsigned int M,
               double* y,
               double* w);

void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const double* maxvals,
               const int* maxes,
               const unsigned int npeaks);
