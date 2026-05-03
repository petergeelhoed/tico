#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "mysync.h"

// Magic number constants
#define TESTSYNCWRITE_ARRAY_SIZE 8
#define TESTSYNCWRITE_LOOP_COUNT 10
int main(void)
{
    struct myarr* testArray = makemyarr(TESTSYNCWRITE_ARRAY_SIZE);
    for (int index = 0; index < TESTSYNCWRITE_LOOP_COUNT; index++)
    {
        testArray->arr[0] = index;
        syncwrite(testArray->arr, TESTSYNCWRITE_ARRAY_SIZE, "testsync");
    }
    freemyarr(testArray);
    wait();
    exit(0);
}
