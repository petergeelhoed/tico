#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "mylib.h"
#include "mysync.h"

#define size 8
int main(void)
{
    struct myarr* test = makemyarr(size);
    for (int i = 0; i < 10; i++)
    {
        test->arr[0] = i;
        syncwrite(test->arr, size, "testsync");
    }
    freemyarr(test);
    wait();
    exit(0);
}
