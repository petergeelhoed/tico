#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "mysync.h"

static pthread_mutex_t count_mutex;

void thread_unlock() { pthread_mutex_unlock(&count_mutex); }

void thread_lock() { pthread_mutex_lock(&count_mutex); }

void writearray(int* arr, unsigned int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    for (unsigned int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %d\n", j, arr[j]);
    }
    fclose(fp);
}

long unsigned int syncappend(int* input, unsigned int NN, FILE* file)
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
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0lu;
    pthread_create(&tid, &attr, threadAppend, info);

    return tid;
}

void* threadAppend(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);
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

    printTOD(file);
    for (unsigned int j = 0; j < mine.NN; j++)
    {
        fprintf(file, "%d\n", mine.array[j]);
    }
    fflush(file);

    free(copyarr);
    pthread_mutex_unlock(&count_mutex);
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

    struct mystruct* info = calloc(1, sizeof *info);

    int* copyarr = calloc(NN, sizeof(int));
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    strcpy(info->file, file);
    info->NN = NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0lu;
    pthread_create(&tid, &attr, threadWrite, info);
}

void* threadWrite(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);
    struct mystruct
    {
        int* array;
        char file[20];
        unsigned int NN;
    } mine = *(struct mystruct*)inStruct;

    int* arrptr = mine.array;
    int* copyarr = calloc(mine.NN, sizeof(int));
    memcpy(copyarr, arrptr, mine.NN * sizeof(int));

    writearray(mine.array, mine.NN, mine.file);

    free(copyarr);
    free(arrptr);
    free(inStruct);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}

void printTOD(FILE* out)
{
    if (out == 0)
        return;
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);

    struct tm* today = localtime(&tv.tv_sec);
    fprintf(out,
            "# %04d-%02d-%02dT%02d:%02d:%02d.%ld %lu.%lu\n",
            today->tm_year + 1900,
            today->tm_mon + 1,
            today->tm_mday,
            today->tm_hour,
            today->tm_min,
            today->tm_sec,
            tv.tv_usec,
            tv.tv_sec,
            tv.tv_usec);
}

void syncappendDouble(double* input, unsigned int NN, FILE* file)
{
    if (file == 0)
        return;
    struct mystruct
    {
        double* array;
        FILE* file;
        unsigned int NN;
    };

    struct mystruct* info = malloc(sizeof *info);

    double* copyarr = malloc(NN * sizeof(double));
    memcpy(copyarr, input, NN * sizeof(double));
    info->array = copyarr;
    info->file = file;
    info->NN = NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0lu;
    pthread_create(&tid, &attr, threadAppendDouble, info);
}

void* threadAppendDouble(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);
    struct mystruct
    {
        double* array;
        FILE* file;
        unsigned int NN;
    } mine = *(struct mystruct*)inStruct;

    double* arrptr = mine.array;
    FILE* file = mine.file;
    double* copyarr = malloc(mine.NN * sizeof(double));
    memcpy(copyarr, arrptr, mine.NN * sizeof(double));
    mine.array = copyarr;
    free(arrptr);
    free(inStruct);

    printTOD(file);
    for (unsigned int j = 0; j < mine.NN; j++)
    {
        fprintf(file, "%f\n", mine.array[j]);
    }
    fflush(file);

    free(copyarr);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}
