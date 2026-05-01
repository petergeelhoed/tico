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
    struct myarr reference = {calloc(ArrayLength, sizeof(int)), 0, ArrayLength};
    struct myarr refdouble = {0,
                              calloc(ArrayLength, sizeof(double)),
                              ArrayLength};
    for (int i = 0; i < (int)ArrayLength; ++i)
    {
        reference.arr[i] = i;
        refdouble.arrd[i] = (double)i * i;
    }
    FILE* file = fopen("out", "a");
    syncAppendMyarr(&reference, file);

    syncAppendMyarr(&refdouble, file);

    waitClose(file);
    return 0;
}
