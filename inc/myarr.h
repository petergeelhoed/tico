#pragma once
#include <stdlib.h>

struct myarr
{
    int* arr;
    double* arrd;
    unsigned int ArrayLength;
};

void freemyarr(struct myarr* arr);
struct myarr* makemyarr(unsigned int ArrayLength);
struct myarr* makemyarrd(unsigned int ArrayLength);
