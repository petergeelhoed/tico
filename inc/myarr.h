#pragma once
#include <stdlib.h>

struct myarr
{
    int* arr;
    double* arrd;
    size_t ArrayLength;
};

void freemyarr(struct myarr* arr);
struct myarr* makemyarr(size_t ArrayLength);
struct myarr* makemyarrd(size_t ArrayLength);
