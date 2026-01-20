#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "mylib.h"
#include "mysync.h"

#define size 8
int main(void)
{
    FILE* filePtr = fopen("testsync", "w");
    struct myarr* test = makemyarr(size);
    printf("%d %d %d \n", test->arr[0], test->arr[1], test->arr[2]);
    for (int i = 0; i < 1000; i++)
    {
        test->arr[0] = i;
        syncAppendMyarr(test, filePtr);
    }
    freemyarr(test);
    wait_close(filePtr);
    exit(0);
}
