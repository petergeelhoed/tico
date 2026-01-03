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
struct myarr* makemyarr(unsigned int NN);
struct myarr* makemyarrd(unsigned int NN);
#endif
