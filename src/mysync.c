
// mysync.c (refactored)

#include <errno.h>
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

/* -----------------------------------------------------------------------------
 * Synchronization primitives
 * -------------------------------------------------------------------------- */

/* Serialize writes to FILE* so output from multiple threads doesn't interleave
 */
static pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;

/* Protect the worker count and provide a way to wait efficiently for "drain" */
static pthread_mutex_t ctr_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t ctr_zero = PTHREAD_COND_INITIALIZER;

static int count = 0;

static void thread_ctr_lock(void) { (void)pthread_mutex_lock(&ctr_mutex); }
static void thread_ctr_unlock(void) { (void)pthread_mutex_unlock(&ctr_mutex); }
static void thread_lock(void) { (void)pthread_mutex_lock(&io_mutex); }
static void thread_unlock(void) { (void)pthread_mutex_unlock(&io_mutex); }

/* -----------------------------------------------------------------------------
 * Task payload passed to worker thread
 * Define once at file scope to avoid UB from multiple anonymous struct tags.
 * -------------------------------------------------------------------------- */
struct append_task
{
    struct myarr* array; /* deep copy of input->arr / input->arrd */
    FILE* file; /* not owned by worker; must stay open until worker finishes */
};

/* -----------------------------------------------------------------------------
 * Blocking wait for all outstanding workers to finish (no busy-wait)
 * -------------------------------------------------------------------------- */
void wait(void) /* keep the original name/signature for compatibility */
{
    thread_ctr_lock();
    while (count > 0)
    {
        (void)pthread_cond_wait(&ctr_zero, &ctr_mutex);
    }
    thread_ctr_unlock();
}

/* -----------------------------------------------------------------------------
 * Utility: write an array to a new file (format-specifier fixed for 'j')
 * -------------------------------------------------------------------------- */
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
        (void)fprintf(
            filePtr, "%u %d\n", j, arr[j]); /* %u for unsigned int j */
    }

    if (fclose(filePtr))
    {
        perror("Error closing file");
    }
}

/* -----------------------------------------------------------------------------
 * Worker thread: serialize writes, flush, free payload, and signal drain
 * -------------------------------------------------------------------------- */
void* threadAppendMyarr(void* inStruct)
{
    struct append_task* mine = (struct append_task*)inStruct;

    /* Serialize all writes to the FILE* to prevent interleaving */
    thread_lock();

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

    thread_unlock();

    /* Free deep copies and task struct */
    free(mine->array->arr);
    free(mine->array->arrd);
    free(mine->array);
    free(mine);

    /* Decrement worker count and signal if we reached zero */
    thread_ctr_lock();
    count--;
    if (count == 0)
    {
        (void)pthread_cond_broadcast(&ctr_zero);
    }
    thread_ctr_unlock();

    return NULL; /* no need to call pthread_exit(NULL) explicitly */
}

/* -----------------------------------------------------------------------------
 * Launch a detached worker that appends 'input' to 'file' safely.
 * Returns: 0 on failure; non-zero token on success (casts pthread_t).
 *
 * NOTE: Casting pthread_t to 'unsigned long' is not fully portable.
 * If you rely on the returned value, consider changing the signature to return
 * 'pthread_t' or simply 'int' for success/failure.
 * -------------------------------------------------------------------------- */
long unsigned int syncAppendMyarr(struct myarr* input, FILE* file)
{
    if (input == NULL || file == NULL)
    {
        errno = EINVAL;
        return 0;
    }

    /* Allocate and deep-copy the payload for the worker thread */
    struct append_task* info = malloc(sizeof(*info));
    if (info == NULL)
    {
        perror("Error allocating task");
        return 0;
    }

    struct myarr* local = malloc(sizeof(*local));
    if (local == NULL)
    {
        perror("Error allocating myarr");
        free(info);
        return 0;
    }
    local->arr = NULL;
    local->arrd = NULL;
    local->ArrayLength = input->ArrayLength;

    /* Copy int array if present and length > 0 */
    if (input->arr != NULL && input->ArrayLength > 0)
    {
        local->arr = malloc(input->ArrayLength * sizeof(int));
        if (local->arr == NULL)
        {
            perror("Error allocating arr");
            free(local);
            free(info);
            return 0;
        }
        memcpy(local->arr, input->arr, input->ArrayLength * sizeof(int));
    }

    /* Copy double array if present and length > 0 */
    if (input->arrd != NULL && input->ArrayLength > 0)
    {
        local->arrd = malloc(input->ArrayLength * sizeof(double));
        if (local->arrd == NULL)
        {
            perror("Error allocating arrd");
            free(local->arr); /* OK if NULL */
            free(local);
            free(info);
            return 0;
        }
        memcpy(local->arrd, input->arrd, input->ArrayLength * sizeof(double));
    }

    info->array = local;
    info->file = file;

    /* Bump worker count before creating the thread */
    thread_ctr_lock();
    count++;
    thread_ctr_unlock();

    /* Create a detached thread */
    pthread_attr_t attr;
    if (pthread_attr_init(&attr) != 0)
    {
        perror("pthread_attr_init failed");

        /* Roll back count */
        thread_ctr_lock();
        count--;
        if (count == 0)
        {
            pthread_cond_broadcast(&ctr_zero);
        }
        thread_ctr_unlock();

        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return 0;
    }
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED) != 0)
    {
        perror("pthread_attr_setdetachstate failed");

        (void)pthread_attr_destroy(&attr); /* clean up attributes */

        /* Roll back count */
        thread_ctr_lock();
        count--;
        if (count == 0)
        {
            pthread_cond_broadcast(&ctr_zero);
        }
        thread_ctr_unlock();

        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return 0;
    }

    pthread_t tid;
    const int succes = pthread_create(&tid, &attr, threadAppendMyarr, info);

    /* Attributes are no longer needed after pthread_create (success or failure)
     */
    (void)pthread_attr_destroy(&attr);

    if (succes != 0)
    {
        perror("Error creating thread");

        /* Roll back count */
        thread_ctr_lock();
        count--;
        if (count == 0)
        {
            pthread_cond_broadcast(&ctr_zero);
        }
        thread_ctr_unlock();

        /* Free payload */
        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return 0;
    }

    return (long unsigned int)tid;
}

void* threadAppend(void* inStruct)
{
    thread_ctr_lock();

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
        thread_ctr_unlock();
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
    thread_ctr_unlock();
    pthread_exit(NULL);
}

void syncwrite(int* input, unsigned int ArrayLength, const char* file)
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
        perror("Error allocating info memory");
        return;
    }

    int* copyarr = calloc(ArrayLength, sizeof(int));
    if (copyarr == NULL)
    {
        perror("Error allocating copyarr memory");
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
    thread_ctr_lock();

    struct mystruct
    {
        int* array;
        char file[FILE_NAME_LENGTH];
        unsigned int ArrayLength;
    }* mine = inStruct;

    writearray(mine->array, mine->ArrayLength, mine->file);

    free(mine->array);
    free(mine);
    thread_ctr_unlock();
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
