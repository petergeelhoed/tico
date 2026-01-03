#include <stdlib.h>
#ifndef MYARR
#define MYARR
struct myarr
{
    int* arr;
    double* arrd;
    unsigned int NN;
};

void freemyarr(struct myarr* arr);
#endif
