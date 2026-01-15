#include "myarr.h"
#include <stdio.h>
void freemyarr(struct myarr* arr)
{
    free(arr->arr);
    free(arr->arrd);
    arr->ArrayLength = 0;
    arr->arr = NULL;
    arr->arrd = NULL;
    free(arr);
}

struct myarr* makemyarrd(unsigned int ArrayLength)
{
    struct myarr* ret = (struct myarr*)calloc(1, sizeof(struct myarr));
    if (ret == NULL)
    {
        return NULL;
    }
    ret->arrd = calloc(ArrayLength, sizeof(double));
    if (ret->arrd == NULL)
    {
        free(ret);
        return NULL;
    }
    ret->ArrayLength = ArrayLength;
    return ret;
}
struct myarr* makemyarr(unsigned int ArrayLength)
{
    struct myarr* ret = (struct myarr*)calloc(1, sizeof(struct myarr));
    if (ret == NULL)
    {
        return NULL;
    }
    ret->arr = calloc(ArrayLength, sizeof(int));
    if (ret->arr == NULL)
    {
        free(ret);
        return NULL;
    }
    ret->ArrayLength = ArrayLength;
    return ret;
}
