#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> // IWYU pragma: keep
#include <unistd.h>

#include "myarr.h"
#include "mydefs.h"
#include "mysync.h"

// NOLINTNEXTLINE(misc-include-cleaner)
static pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
volatile int count = 0;

void wait(void)
{
    while (count > 0)
    {
        const int WAIT = 10000;
        usleep(WAIT);
    }
}

void thread_unlock(void) { pthread_mutex_unlock(&count_mutex); }

void thread_lock(void) { pthread_mutex_lock(&count_mutex); }

void writearray(int* arr, unsigned int ArrayLength, const char* file)
{
    FILE* filePtr = fopen(file, "w");
    if (filePtr == NULL)
    {
        perror("Error opening file");
        return;
    }
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        (void)fprintf(filePtr, "%d %d\n", j, arr[j]);
    }
    if (fclose(filePtr))
    {
        perror("Error closing file");
    }
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
        free(info);
        count--;
        perror("Error allocating memory");
        return 0;
    }
    info->array->arr = NULL;
    info->array->arrd = NULL;

    if (input->arr != NULL)
    {
        info->array->arr = malloc(input->ArrayLength * sizeof(int));
        if (info->array->arr == NULL)
        {
            free(info);
            count--;
            perror("Error allocating memory");
            return 0;
        }

        memcpy(info->array->arr, input->arr, input->ArrayLength * sizeof(int));
    }
    if (input->arrd != NULL)
    {
        info->array->arrd = malloc(input->ArrayLength * sizeof(double));
        if (info->array->arrd == NULL)
        {
            free(info);
            count--;
            perror("Error allocating memory");
            return 0;
        }

        memcpy(info->array->arrd,
               input->arrd,
               input->ArrayLength * sizeof(double));
    }

    info->file = file;
    info->array->ArrayLength = input->ArrayLength;

    pthread_attr_t attr; // NOLINT(misc-include-cleaner)
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_t tid = 0; // NOLINT(misc-include-cleaner)
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
        for (unsigned int j = 0; j < mine->array->ArrayLength; j++)
        {
            (void)fprintf(mine->file, "%d\n", mine->array->arr[j]);
        }
    }

    if (mine->array->arrd != NULL)
    {
        for (unsigned int j = 0; j < mine->array->ArrayLength; j++)
        {
            (void)fprintf(mine->file, "%f\n", mine->array->arrd[j]);
        }
    }

    (void)fflush(mine->file);

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
        unsigned int ArrayLength;
    }* mine = inStruct;

    int* copyarr = malloc(mine->ArrayLength * sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        pthread_mutex_unlock(&count_mutex);
        pthread_exit(NULL);
    }
    memcpy(copyarr, mine->array, mine->ArrayLength * sizeof(int));
    free(mine->array);
    mine->array = copyarr;

    printTOD(mine->file);
    for (unsigned int j = 0; j < mine->ArrayLength; j++)
    {
        (void)fprintf(mine->file, "%d\n", mine->array[j]);
    }
    (void)fflush(mine->file);

    free(mine->array);
    free(mine);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}

void syncwrite(int* input, unsigned int ArrayLength, char* file)
{
    struct mystruct
    {
        int* array;
        char file[FILE_NAME_LENGTH];
        unsigned int ArrayLength;
    };

    struct mystruct* info = calloc(1, sizeof *info);
    if (info == NULL)
    {
        perror("Error allocating memory");
        return;
    }

    int* copyarr = calloc(ArrayLength, sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating memory");
        free(info);
        return;
    }
    memcpy(copyarr, input, ArrayLength * sizeof(int));
    info->array = copyarr;
    strncpy(info->file, file, FILE_NAME_LENGTH - 1);
    info->file[FILE_NAME_LENGTH - 1] = '\0';
    info->ArrayLength = ArrayLength;

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
        unsigned int ArrayLength;
    }* mine = inStruct;

    writearray(mine->array, mine->ArrayLength, mine->file);

    free(mine->array);
    free(mine);
    pthread_mutex_unlock(&count_mutex);
    pthread_exit(NULL);
}

void printTOD(FILE* out)
{
    if (out == NULL)
    {
        return;
    }

    struct timeval time; // NOLINT(misc-include-cleaner)
    gettimeofday(&time, NULL);

    struct tm* today = localtime(&time.tv_sec); // NOLINT(misc-include-cleaner)
    if (today == NULL)
    {
        perror("Error getting local time");
        return;
    }

    const int nineteenhundred = 1900;
    (void)fprintf(out,
                  "# %04d-%02d-%02dT%02d:%02d:%02d.%06ld %lu.%06lu\n",
                  today->tm_year + nineteenhundred,
                  today->tm_mon + 1,
                  today->tm_mday,
                  today->tm_hour,
                  today->tm_min,
                  today->tm_sec,
                  time.tv_usec,
                  time.tv_sec,
                  time.tv_usec);
}
