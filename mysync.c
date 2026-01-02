#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "myarr.h"
#include "mysync.h"

#define FILE_NAME_LENGTH 256

static pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
volatile int count = 0;

void wait()
{
    while (count > 0)
    {
        usleep(10000);
    }
}

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

long unsigned int syncAppendMyarr(struct myarr* input, FILE* file)
{
    count++;
    struct mystruct
    {
        struct myarr* array;
        FILE* file;
    };

    struct mystruct* info = malloc(sizeof(struct mystruct));
    if (info == NULL)
    {
        count--;
        perror("Error allocating memory");
        return 0;
    }

    info->array = malloc(sizeof(struct myarr));
    if (info->array == NULL)
    {
        count--;
        perror("Error allocating memory");
        return 0;
    }
    info->array->arr = NULL;
    info->array->arrd = NULL;

    if (input->arr != NULL)
    {
        info->array->arr = malloc(input->NN * sizeof(int));
        if (info->array->arr == NULL)
        {
            count--;
            perror("Error allocating memory");
            return 0;
        }

        memcpy(info->array->arr, input->arr, input->NN * sizeof(int));
    }
    if (input->arrd != NULL)
    {
        info->array->arrd = malloc(input->NN * sizeof(double));
        if (info->array->arrd == NULL)
        {
            count--;
            perror("Error allocating memory");
            return 0;
        }

        memcpy(info->array->arrd, input->arrd, input->NN * sizeof(double));
    }

    info->file = file;
    info->array->NN = input->NN;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0;
    if (pthread_create(&tid, &attr, threadAppendMyarr, info) != 0)
    {
        perror("Error creating thread");
        free(info->array->arr);
        free(info->array->arrd);
        free(info);
        return 0;
    }
    return tid;
}

void* threadAppendMyarr(void* inStruct)
{
    pthread_mutex_lock(&count_mutex);

    struct mystruct
    {
        struct myarr* array;
        FILE* file;
    }* mine = inStruct;

    printTOD(mine->file);
    if (mine->array->arr != NULL)
    {
        for (unsigned int j = 0; j < mine->array->NN; j++)
        {
            fprintf(mine->file, "%d\n", mine->array->arr[j]);
        }
    }

    if (mine->array->arrd != NULL)
    {
        for (unsigned int j = 0; j < mine->array->NN; j++)
        {
            fprintf(mine->file, "%f\n", mine->array->arrd[j]);
        }
    }

    fflush(mine->file);

    free(mine->array->arr);
    free(mine->array->arrd);
    free(mine->array);
    free(mine);
    count--;
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
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
