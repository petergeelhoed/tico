#include <errno.h>
#include <limits.h>
#include <pthread.h>
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

/** @brief Tracks pending writes and close state for one output file. */
struct file_state
{
    /** File associated with this state. */
    FILE* file;
    /** Number of append operations still pending for the file. */
    unsigned int pending;
    /** Nonzero once closing has been requested. */
    int closing;
    /** Next file state in the global linked list. */
    struct file_state* next;
};

static struct file_state* file_states = NULL;

static void report_errno_error(const char* operation)
{
    const int error_code = errno;
    (void)fprintf(stderr, "%s: %s\n", operation, strerror(error_code));
}

static void report_pthread_error(const char* operation, int error_code)
{
    if (error_code != 0)
    {
        (void)fprintf(stderr, "%s: %s\n", operation, strerror(error_code));
    }
}

/**
 * @brief Lock the counter mutex.
 */
static void thread_ctr_lock(void)
{
    report_pthread_error("pthread_mutex_lock", pthread_mutex_lock(&ctr_mutex));
}
/**
 * @brief Unlock the counter mutex.
 */
static void thread_ctr_unlock(void)
{
    report_pthread_error("pthread_mutex_unlock",
                         pthread_mutex_unlock(&ctr_mutex));
}
/**
 * @brief Lock the IO mutex.
 */
static void thread_lock(void)
{
    report_pthread_error("pthread_mutex_lock", pthread_mutex_lock(&io_mutex));
}
/**
 * @brief Unlock the IO mutex.
 */
static void thread_unlock(void)
{
    report_pthread_error("pthread_mutex_unlock",
                         pthread_mutex_unlock(&io_mutex));
}

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

/**
 * @brief Remove the file state entry for a given file.
 *
 * @param file The file pointer whose state should be removed.
 */
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
static void printTODUnlockedTime(FILE* out, struct timeval* time)
{
    if (out == NULL || time == NULL)
    {
        return;
    }
    struct tm now;
    if (localtime_r(&time->tv_sec, &now) == NULL)
    {
        report_errno_error("localtime_r");
        return;
    }

    const int nineteenhundred = 1900;
    if (fprintf(out,
                "# %04d-%02d-%02dT%02d:%02d:%02d.%06ld %lld.%06ld\n",
                now.tm_year + nineteenhundred,
                now.tm_mon + 1,
                now.tm_mday,
                now.tm_hour,
                now.tm_min,
                now.tm_sec,
                (long)time->tv_usec,
                (long long)time->tv_sec,
                (long)time->tv_usec) < 0)
    {
        perror("printTODUnlockedTime: fprintf");
    }
}

/**
 * @brief Print the current date and time to the given output file.
 *
 * @param out Output file pointer.
 */
static void printTODUnlocked(FILE* out)
{
    if (out == NULL)
    {
        return;
    }

    struct timeval time;
    if (gettimeofday(&time, NULL) != 0)
    {
        report_errno_error("printTODUnlocked: gettimeofday");
        return;
    }
    printTODUnlockedTime(out, &time);
}

/**
 * @brief Decrement the global count and signal waiting threads if needed.
 */
static void decr_count(void)
{
    thread_ctr_lock();
    count--;
    report_pthread_error("pthread_cond_broadcast",
                         pthread_cond_broadcast(&ctr_zero));
    thread_ctr_unlock();
}

/**
 * @brief Decrement the append count for a file and signal waiting threads if
 * needed.
 *
 * @param file The file pointer whose append count should be decremented.
 */
static void decr_append_count(FILE* file)
{
    thread_ctr_lock();
    count--;

    struct file_state* state = find_file_state(file);
    if (state != NULL && state->pending > 0)
    {
        state->pending--;
    }

    report_pthread_error("pthread_cond_broadcast",
                         pthread_cond_broadcast(&ctr_zero));
    thread_ctr_unlock();
}
/** @brief Data passed to a worker that appends an array to a file. */
struct append_task
{
    /** Deep copy of the input array. */
    struct myarr* array;
    /** Output file, which must remain open until the worker finishes. */
    FILE* file;
    /** Timestamp written before the array data. */
    struct timeval time;
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
        perror("waitClose");
        return;
    }

    state->closing = 1;
    while (state->pending > 0)
    {
        int error_code = pthread_cond_wait(&ctr_zero, &ctr_mutex);
        if (error_code != 0)
        {
            report_pthread_error("pthread_cond_wait", error_code);
            thread_ctr_unlock();
            return;
        }
    }

    remove_file_state(file);
    thread_ctr_unlock();

    if (fclose(file) != 0)
    {
        perror("waitClose: fclose");
    }
}

void wait(void)
{
    thread_ctr_lock();
    while (count > 0)
    {
        int error_code = pthread_cond_wait(&ctr_zero, &ctr_mutex);
        if (error_code != 0)
        {
            report_pthread_error("pthread_cond_wait", error_code);
            thread_ctr_unlock();
            return;
        }
    }
    thread_ctr_unlock();
}

void* threadAppendMyarr(void* inStruct)
{
    struct append_task* mine = (struct append_task*)inStruct;
    FILE* file = mine->file;

    /* Serialize all writes to the FILE* to prevent interleaving */
    thread_lock();

    printTODUnlockedTime(mine->file, &mine->time);

    if (mine->array->arr != NULL)
    {
        for (unsigned int j = 0; j < mine->array->ArrayLength; j++)
        {
            if (fprintf(mine->file, "%d\n", mine->array->arr[j]) < 0)
            {
                perror("threadAppendMyarr: fprintf");
                break;
            }
        }
    }
    if (mine->array->arrd != NULL)
    {
        for (unsigned int j = 0; j < mine->array->ArrayLength; j++)
        {
            if (fprintf(mine->file, "%f\n", mine->array->arrd[j]) < 0)
            {
                perror("threadAppendMyarr: fprintf");
                break;
            }
        }
    }

    if (fflush(mine->file) != 0)
    {
        perror("threadAppendMyarr: fflush");
    }

    thread_unlock();

    /* Free deep copies and task struct */
    freemyarr(mine->array);
    free(mine);

    decr_append_count(file);
    return NULL;
}

int syncAppendMyarr(struct myarr* input, FILE* file)
{
    if (input == NULL || file == NULL)
    {
        errno = EINVAL;
        perror("syncAppendMyarr");
        return -1;
    }

    struct timeval time;
    if (gettimeofday(&time, NULL) != 0)
    {
        report_errno_error("syncAppendMyarr: gettimeofday");
        return -1;
    }

    struct append_task* info = malloc(sizeof(*info));
    if (info == NULL)
    {
        perror("Error allocating task");
        return -1;
    }

    struct myarr* local = input->arr != NULL ? makemyarr(input->ArrayLength)
                                             : makemyarrd(input->ArrayLength);
    if (local == NULL)
    {
        perror("Error allocating myarr");
        free(info);
        return -1;
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
    memcpy(&info->time, &time, sizeof(struct timeval));

    thread_ctr_lock();
    struct file_state* state = get_or_create_file_state(file);
    if (state == NULL)
    {
        thread_ctr_unlock();
        perror("Error tracking file state");
        freemyarr(local);
        free(info);
        return -1;
    }
    if (state->closing)
    {
        thread_ctr_unlock();
        errno = EPIPE;
        perror("syncAppendMyarr: file is closing");
        freemyarr(local);
        free(info);
        return -1;
    }
    state->pending++;
    count++;
    thread_ctr_unlock();

    pthread_attr_t attr;
    int error_code = pthread_attr_init(&attr);
    if (error_code != 0)
    {
        report_pthread_error("pthread_attr_init", error_code);

        decr_append_count(file);

        freemyarr(local);
        free(info);
        return -1;
    }
    error_code = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    if (error_code != 0)
    {
        report_pthread_error("pthread_attr_setdetachstate", error_code);

        report_pthread_error("pthread_attr_destroy",
                             pthread_attr_destroy(&attr));

        decr_append_count(file);

        freemyarr(local);
        free(info);
        return -1;
    }

    pthread_t tid;
    const int succes = pthread_create(&tid, &attr, threadAppendMyarr, info);
    report_pthread_error("pthread_attr_destroy", pthread_attr_destroy(&attr));

    if (succes != 0)
    {
        report_pthread_error("pthread_create", succes);

        decr_append_count(file);

        /* Free payload */
        freemyarr(local);
        free(info);
        return -1;
    }

    return 0;
}

/** @brief Data passed to a worker that writes integer values to a file. */
struct write_task
{
    /** Owned copy of the values to write. */
    int* array;
    /** Number of values in array. */
    size_t len;
    /** Output file path. */
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
    if (file != NULL)
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
    else
    {
        perror("threadWrite: fopen");
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
        perror("syncwrite");
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

    thread_ctr_lock();
    count++;
    thread_ctr_unlock();

    pthread_attr_t attr;
    int error_code = pthread_attr_init(&attr);
    if (error_code != 0)
    {
        report_pthread_error("pthread_attr_init", error_code);
        decr_count();
        free(copyarr);
        free(task);
        return -1;
    }
    error_code = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    if (error_code != 0)
    {
        report_pthread_error("pthread_attr_setdetachstate", error_code);
        report_pthread_error("pthread_attr_destroy",
                             pthread_attr_destroy(&attr));
        decr_count();
        free(copyarr);
        free(task);
        return -1;
    }

    pthread_t tid;
    int success = pthread_create(&tid, &attr, threadWrite, task);
    report_pthread_error("pthread_attr_destroy", pthread_attr_destroy(&attr));

    if (success != 0)
    {
        report_pthread_error("pthread_create", success);
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
