#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "mysync.h"

#define FILE_NAME_LENGTH 256

static pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;

void thread_unlock() { pthread_mutex_unlock(&count_mutex); }

void thread_lock() { pthread_mutex_lock(&count_mutex); }

void writearray(int* arr, unsigned int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    if (fp == NULL)
    {
        perror("Error opening file");
        return;
    }
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
    if (info == NULL)
    {
        perror("Error allocating memory");
        return 0;
    }

    int* copyarr = malloc(NN * sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        free(info);
        return 0;
    }
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    info->file = file;
    info->NN = NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0;
    if (pthread_create(&tid, &attr, threadAppend, info) != 0)
    {
        perror("Error creating thread");
        free(copyarr);
        free(info);
        return 0;
    }

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
    }* mine = inStruct;

    int* copyarr = malloc(mine->NN * sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        pthread_mutex_unlock(&count_mutex);
        pthread_exit(NULL);
    }
    memcpy(copyarr, mine->array, mine->NN * sizeof(int));
    free(mine->array);
    mine->array = copyarr;

    printTOD(mine->file);
    for (unsigned int j = 0; j < mine->NN; j++)
    {
        fprintf(mine->file, "%d\n", mine->array[j]);
    }
    fflush(mine->file);

    free(mine->array);
    free(mine);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}

void syncwrite(int* input, unsigned int NN, char* file)
{
    struct mystruct
    {
        int* array;
        char file[FILE_NAME_LENGTH];
        unsigned int NN;
    };

    struct mystruct* info = calloc(1, sizeof *info);
    if (info == NULL)
    {
        perror("Error allocating memory");
        return;
    }

    int* copyarr = calloc(NN, sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        free(info);
        return;
    }
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    strncpy(info->file, file, FILE_NAME_LENGTH - 1);
    info->file[FILE_NAME_LENGTH - 1] = '\0';
    info->NN = NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0;
    if (pthread_create(&tid, &attr, threadWrite, info) != 0)
    {
        perror("Error creating thread");
        free(copyarr);
        free(info);
    }
}

void* threadWrite(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);

    struct mystruct
    {
        int* array;
        char file[FILE_NAME_LENGTH];
        unsigned int NN;
    }* mine = inStruct;

    writearray(mine->array, mine->NN, mine->file);

    free(mine->array);
    free(mine);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}

void printTOD(FILE* out)
{
    if (out == NULL)
        return;

    struct timeval tv;
    gettimeofday(&tv, NULL);

    struct tm* today = localtime(&tv.tv_sec);
    if (today == NULL)
    {
        perror("Error getting local time");
        return;
    }

    fprintf(out,
            "# %04d-%02d-%02dT%02d:%02d:%02d.%06ld %lu.%06lu\n",
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
    if (file == NULL)
        return;

    struct mystruct
    {
        double* array;
        FILE* file;
        unsigned int NN;
    };

    struct mystruct* info = malloc(sizeof *info);
    if (info == NULL)
    {
        perror("Error allocating memory");
        return;
    }

    double* copyarr = malloc(NN * sizeof(double));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        free(info);
        return;
    }
    memcpy(copyarr, input, NN * sizeof(double));
    info->array = copyarr;
    info->file = file;
    info->NN = NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0;
    if (pthread_create(&tid, &attr, threadAppendDouble, info) != 0)
    {
        perror("Error creating thread");
        free(copyarr);
        free(info);
    }
}

void* threadAppendDouble(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);

    struct mystruct
    {
        double* array;
        FILE* file;
        unsigned int NN;
    }* mine = inStruct;

    double* copyarr = malloc(mine->NN * sizeof(double));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        pthread_mutex_unlock(&count_mutex);
        pthread_exit(NULL);
    }
    memcpy(copyarr, mine->array, mine->NN * sizeof(double));
    free(mine->array);
    mine->array = copyarr;

    printTOD(mine->file);
    for (unsigned int j = 0; j < mine->NN; j++)
    {
        fprintf(mine->file, "%f\n", mine->array[j]);
    }
    fflush(mine->file);

    free(mine->array);
    free(mine);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}
