#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "mylib.h"
#include "mysync.h"

// Magic number constants
#define TESTSYNC_SIZE 8
#define TESTSYNC_LOOP_COUNT 1000
int main(void)
{
    FILE* filePtr = fopen("testsync", "w");
    struct myarr* test = makemyarr(TESTSYNC_SIZE);
    printf("%d %d %d \n", test->arr[0], test->arr[1], test->arr[2]);
    for (int i = 0; i < TESTSYNC_LOOP_COUNT; i++)
    {
        test->arr[0] = i;
        syncAppendMyarr(test, filePtr);
    }
    freemyarr(test);
    waitClose(filePtr);
    exit(0);
}
