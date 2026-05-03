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

static pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t ctr_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t ctr_zero = PTHREAD_COND_INITIALIZER;

static int count = 0;

struct file_state
{
    FILE* file;
    unsigned int pending;
    int closing;
    struct file_state* next;
};

static struct file_state* file_states = NULL;

static void thread_ctr_lock(void) { (void)pthread_mutex_lock(&ctr_mutex); }
static void thread_ctr_unlock(void) { (void)pthread_mutex_unlock(&ctr_mutex); }
static void thread_lock(void) { (void)pthread_mutex_lock(&io_mutex); }
static void thread_unlock(void) { (void)pthread_mutex_unlock(&io_mutex); }

static struct file_state* find_file_state(FILE* file)
{
    for (struct file_state* it = file_states; it != NULL; it = it->next)
    {
        if (it->file == file)
        {
            return it;
        }
    }
    return NULL;
}

static struct file_state* get_or_create_file_state(FILE* file)
{
    struct file_state* state = find_file_state(file);
    if (state != NULL)
    {
        return state;
    }

    state = calloc(1, sizeof(*state));
    if (state == NULL)
    {
        return NULL;
    }
    state->file = file;
    state->next = file_states;
    file_states = state;
    return state;
}

static void remove_file_state(FILE* file)
{
    struct file_state* prev = NULL;
    struct file_state* iter = file_states;
    while (iter != NULL)
    {
        if (iter->file == file)
        {
            if (prev == NULL)
            {
                file_states = iter->next;
            }
            else
            {
                prev->next = iter->next;
            }
            free(iter);
            return;
        }
        prev = iter;
        iter = iter->next;
    }
}

static void printTODUnlocked(FILE* out)
{
    if (out == NULL)
    {
        return;
    }

    struct timeval time;
    gettimeofday(&time, NULL);

    struct tm today;
    if (localtime_r(&time.tv_sec, &today) == NULL)
    {
        perror("Error getting local time");
        return;
    }

    const int nineteenhundred = 1900;
    (void)fprintf(out,
                  "# %04d-%02d-%02dT%02d:%02d:%02d.%06ld %lld.%06ld\n",
                  today.tm_year + nineteenhundred,
                  today.tm_mon + 1,
                  today.tm_mday,
                  today.tm_hour,
                  today.tm_min,
                  today.tm_sec,
                  (long)time.tv_usec,
                  (long long)time.tv_sec,
                  (long)time.tv_usec);
}

static void decr_count(void)
{
    thread_ctr_lock();
    count--;
    pthread_cond_broadcast(&ctr_zero);
    thread_ctr_unlock();
}

static void decr_append_count(FILE* file)
{
    thread_ctr_lock();
    count--;

    struct file_state* state = find_file_state(file);
    if (state != NULL && state->pending > 0)
    {
        state->pending--;
    }

    pthread_cond_broadcast(&ctr_zero);
    thread_ctr_unlock();
}
struct append_task
{
    struct myarr* array; /* deep copy of input->arr / input->arrd */
    FILE* file; /* not owned by worker; must stay open until worker finishes */
};

void waitClose(FILE* file)
{
    if (file == NULL)
    {
        return;
    }

    thread_ctr_lock();
    struct file_state* state = get_or_create_file_state(file);
    if (state == NULL)
    {
        thread_ctr_unlock();
        errno = ENOMEM;
        return;
    }

    state->closing = 1;
    while (state->pending > 0)
    {
        (void)pthread_cond_wait(&ctr_zero, &ctr_mutex);
    }

    remove_file_state(file);
    thread_ctr_unlock();

    (void)fclose(file);
}

void wait(void)
{
    thread_ctr_lock();
    while (count > 0)
    {
        (void)pthread_cond_wait(&ctr_zero, &ctr_mutex);
    }
    thread_ctr_unlock();
}

void* threadAppendMyarr(void* inStruct)
{
    struct append_task* mine = (struct append_task*)inStruct;
    FILE* file = mine->file;

    /* Serialize all writes to the FILE* to prevent interleaving */
    thread_lock();

    printTODUnlocked(mine->file);

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

    decr_append_count(file);
    return NULL;
}

void syncAppendMyarr(struct myarr* input, FILE* file)
{
    if (input == NULL || file == NULL)
    {
        errno = EINVAL;
        return;
    }

    struct append_task* info = malloc(sizeof(*info));
    if (info == NULL)
    {
        perror("Error allocating task");
        return;
    }

    struct myarr* local = input->arr != NULL ? makemyarr(input->ArrayLength)
                                             : makemyarrd(input->ArrayLength);
    if (local == NULL)
    {
        perror("Error allocating myarr");
        free(info);
        return;
    }
    if (input->arr != NULL && input->ArrayLength > 0)
    {
        memcpy(local->arr, input->arr, input->ArrayLength * sizeof(int));
    }

    if (input->arrd != NULL && input->ArrayLength > 0)
    {
        memcpy(local->arrd, input->arrd, input->ArrayLength * sizeof(double));
    }

    info->array = local;
    info->file = file;

    thread_ctr_lock();
    struct file_state* state = get_or_create_file_state(file);
    if (state == NULL)
    {
        thread_ctr_unlock();
        perror("Error tracking file state");
        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return;
    }
    if (state->closing)
    {
        thread_ctr_unlock();
        errno = EPIPE;
        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return;
    }
    state->pending++;
    count++;
    thread_ctr_unlock();

    pthread_attr_t attr;
    if (pthread_attr_init(&attr) != 0)
    {
        perror("pthread_attr_init failed");

        decr_append_count(file);

        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return;
    }
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED) != 0)
    {
        perror("pthread_attr_setdetachstate failed");

        (void)pthread_attr_destroy(&attr);

        decr_append_count(file);

        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
        return;
    }

    pthread_t tid;
    const int succes = pthread_create(&tid, &attr, threadAppendMyarr, info);
    (void)pthread_attr_destroy(&attr);

    if (succes != 0)
    {
        perror("Error creating thread");

        decr_append_count(file);

        /* Free payload */
        free(local->arrd);
        free(local->arr);
        free(local);
        free(info);
    }
}

struct write_task
{
    int* array;
    size_t len;
    char file[FILE_NAME_LENGTH];
};

static void* threadWrite(void* arg)
{
    struct write_task* task = (struct write_task*)arg;
    if (task == NULL)
    {
        decr_count();
        return NULL;
    }
    FILE* file = fopen(task->file, "w");
    if (file)
    {

        for (unsigned int i = 0; i < task->len; ++i)
        {
            if (fprintf(file, "%d\n", task->array[i]) < 0)
            {
                perror("threadWrite: fprintf failed");
                break;
            }
        }

        if (fflush(file) != 0)
        {
            perror("threadWrite: fflush failed");
        }
        if (fclose(file) != 0)
        {
            perror("threadWrite: fclose failed");
        }
    }
    free(task ? task->array : NULL);
    free(task);
    decr_count();
    return NULL;
}

int syncwrite(int* input, size_t ArrayLength, const char* file)
{
    if (!input || !file)
    {
        errno = EINVAL;
        return -1;
    }

    struct write_task* task = calloc(1, sizeof *task);
    if (!task)
    {
        perror("alloc task");
        return -1;
    }

    int* copyarr = NULL;
    if (ArrayLength > 0)
    {
        copyarr = calloc(ArrayLength, sizeof *copyarr);
        if (!copyarr)
        {
            perror("alloc copyarr");
            free(task);
            return -1;
        }
        memcpy(copyarr, input, ArrayLength * sizeof *copyarr);
    }
    task->array = copyarr;
    task->len = ArrayLength;
    strncpy(task->file, file, FILE_NAME_LENGTH - 1);
    task->file[FILE_NAME_LENGTH - 1] = '\0';

    pthread_mutex_lock(&ctr_mutex);
    count++;
    pthread_mutex_unlock(&ctr_mutex);

    pthread_attr_t attr;
    if (pthread_attr_init(&attr) != 0)
    {
        perror("attr_init");
        decr_count();
        free(copyarr);
        free(task);
        return -1;
    }
    if (pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED) != 0)
    {
        perror("setdetach");
        pthread_attr_destroy(&attr);
        decr_count();
        free(copyarr);
        free(task);
        return -1;
    }

    pthread_t tid;
    int success = pthread_create(&tid, &attr, threadWrite, task);
    pthread_attr_destroy(&attr);

    if (success != 0)
    {
        perror("pthread_create");
        decr_count();
        free(copyarr);
        free(task);
        return -1;
    }
    return 0;
}

void printTOD(FILE* out)
{
    thread_lock();
    printTODUnlocked(out);
    thread_unlock();
}
