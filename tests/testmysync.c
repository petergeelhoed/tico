#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myarr.h"
#include "mysync.h"

int main(void)
{
#define TESTMYSYNC_ARRAY_LENGTH 7
    unsigned int ArrayLength = TESTMYSYNC_ARRAY_LENGTH;
    struct myarr* reference = makemyarr(ArrayLength);
    struct myarr* refdouble = makemyarrd(ArrayLength);
    for (int i = 0; i < (int)ArrayLength; ++i)
    {
        reference->arr[i] = i;
        refdouble->arrd[i] = (double)i * i;
    }
    FILE* file = fopen("out", "a");
    if (file == NULL || syncAppendMyarr(reference, file) != 0 ||
        syncAppendMyarr(refdouble, file) != 0)
    {
        if (file != NULL)
        {
            waitClose(file);
        }
        freemyarr(reference);
        freemyarr(refdouble);
        return EXIT_FAILURE;
    }

    if (syncAppendMyarr(reference, NULL) != -1)
    {
        waitClose(file);
        freemyarr(reference);
        freemyarr(refdouble);
        return EXIT_FAILURE;
    }

    waitClose(file);
    freemyarr(reference);
    freemyarr(refdouble);
    return 0;
}
