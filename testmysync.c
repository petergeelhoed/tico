#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myarr.h"
#include "mysync.h"

int main()
{
    unsigned int NN = 7;
    struct myarr reference = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr refdouble = {0, calloc(NN, sizeof(double)), NN};
    for (unsigned int i = 0; i < NN; ++i)
    {
        reference.arr[i] = i;
        refdouble.arrd[i] = (double)i * i;
    }
    FILE* file = fopen("out", "a");
    syncAppendMyarr(&reference, file);

    syncAppendMyarr(&refdouble, file);

    wait();
    thread_lock();
    fclose(file);
    thread_unlock();
    return 0;
}
