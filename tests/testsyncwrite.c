#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "mylib.h"
#include "mysync.h"

#define kArraySize 8
int main(void)
{
    struct myarr* testArray = makemyarr(kArraySize);
    for (int index = 0; index < 10; index++) // NOLINT(readability-magic-numbers)
    {
        testArray->arr[0] = index;
        syncwrite(testArray->arr, kArraySize, "testsync");
    }
    freemyarr(testArray);
    wait();
    exit(0);
}
