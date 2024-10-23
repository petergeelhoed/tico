#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mysync.h"

void writearray(int* arr, unsigned int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    for (unsigned int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %d\n", j, arr[j]);
    }
    fclose(fp);
}

void syncappend(int* input, unsigned int NN, FILE* file)
{
    struct mystruct
    {
        int* array;
        FILE* file;
        unsigned int NN;
    };

    struct mystruct* info = malloc(sizeof *info);

    int* copyarr = malloc(NN * sizeof(int));
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    info->file = file;
    info->NN = NN;

    pthread_t tid;
    pthread_create(&tid, NULL, threadAppend, info);
    pthread_detach(tid);
}

void* threadAppend(void* inStruct)
{
    struct mystruct
    {
        int* array;
        FILE* file;
        unsigned int NN;
    } mine = *(struct mystruct*)inStruct;

    int* arrptr = mine.array;
    FILE* file = mine.file;
    int* copyarr = malloc(mine.NN * sizeof(int));
    memcpy(copyarr, arrptr, mine.NN * sizeof(int));
    mine.array = copyarr;
    free(arrptr);
    free(inStruct);

    for (unsigned int j = 0; j < mine.NN; j++)
    {
        fprintf(file, "%d\n", mine.array[j]);
    }
    fflush(file);

    free(copyarr);
    pthread_exit(NULL);
}

void syncwrite(int* input, unsigned int NN, char* file)
{
    struct mystruct
    {
        int* array;
        char file[20];
        unsigned int NN;
    };

    struct mystruct* info = malloc(sizeof *info);

    int* copyarr = malloc(NN * sizeof(int));
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    strcpy(info->file, file);
    info->NN = NN;

    pthread_t tid;
    pthread_create(&tid, NULL, threadWrite, info);
    pthread_detach(tid);
}

void* threadWrite(void* inStruct)
{
    struct mystruct
    {
        int* array;
        char file[20];
        unsigned int NN;
    } mine = *(struct mystruct*)inStruct;

    int* arrptr = mine.array;
    int* copyarr = malloc(mine.NN * sizeof(int));
    memcpy(copyarr, arrptr, mine.NN * sizeof(int));
    mine.array = copyarr;
    free(arrptr);
    free(inStruct);

    writearray(mine.array, mine.NN, mine.file);

    pthread_exit(NULL);
}
