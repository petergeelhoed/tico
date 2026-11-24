#include "myarr.h"
#include "printing.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define INIT_ARRSIZE 10
#define BUFFER_SIZE 512

void freemyarr(struct myarr* arr)
{
    free(arr->arr);
    free(arr->arrd);
    arr->capacity = 0;
    arr->ArrayLength = 0;
    arr->arr = NULL;
    arr->arrd = NULL;
    free(arr);
}

struct myarr* makemyarrd(size_t ArrayLength)
{
    struct myarr* ret = (struct myarr*)calloc(1, sizeof(struct myarr));
    if (ret == NULL)
    {
        print("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    ret->arrd = calloc(ArrayLength, sizeof(double));
    if (ret->arrd == NULL)
    {
        free(ret);
        print("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    ret->ArrayLength = ArrayLength;
    ret->capacity = ArrayLength;
    return ret;
}
struct myarr* makemyarr(size_t ArrayLength)
{
    struct myarr* ret = (struct myarr*)calloc(1, sizeof(struct myarr));
    if (ret == NULL)
    {
        print("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    ret->arr = calloc(ArrayLength, sizeof(int));
    if (ret->arr == NULL)
    {
        free(ret);
        print("Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }
    ret->ArrayLength = ArrayLength;
    ret->capacity = ArrayLength;
    return ret;
}

/**
 * @brief Safely reallocates memory for a double array, exits on failure.
 *
 * @param ptr Pointer to the current memory block.
 * @param new_size New size for the array.
 * @return Pointer to the reallocated memory.
 */
static double* safeRealloc(double* ptr, size_t new_size)
{
    double* next = realloc(ptr, new_size * sizeof(double));
    if (!next)
    {
        print("Memory allocation failed\n");
        free(ptr);
        exit(EXIT_FAILURE);
    }
    return next;
}

static int* safeReallocInt(int* ptr, size_t new_size)
{
    int* next = realloc(ptr, new_size * sizeof(int));
    if (!next)
    {
        print("Memory allocation failed\n");
        free(ptr);
        exit(EXIT_FAILURE);
    }
    return next;
}

struct myarr* readinput(void)
{
    struct myarr* arrd = makemyarrd(INIT_ARRSIZE);
    char line[BUFFER_SIZE];
    char* ptr;
    char* endptr;

    while (fgets(line, sizeof(line), stdin))
    {
        ptr = line;
        while (*ptr)
        {
            errno = 0;
            double readDouble = strtod(ptr, &endptr);
            if (ptr == endptr)
            {
                ptr++;
                continue;
            }

            if (arrd->ArrayLength >= arrd->capacity)
            {
                arrd->capacity = (arrd->capacity * 3) / 2;
                arrd->arrd = safeRealloc(arrd->arrd, arrd->capacity);
            }

            arrd->arrd[arrd->ArrayLength++] = readDouble;
            ptr = endptr;
        }
    }

    return arrd;
}

struct myarr* readinputint(void)
{
    struct myarr* arr = makemyarr(INIT_ARRSIZE);
    char line[BUFFER_SIZE];
    char* ptr;
    char* endptr;

    while (fgets(line, sizeof(line), stdin))
    {
        ptr = line;
        while (*ptr)
        {
            errno = 0;
            const int DECIMAL_BASE = 10;
            long readLong = strtol(ptr, &endptr, DECIMAL_BASE);
            if (ptr == endptr)
            {
                ptr++;
                continue;
            }

            if (arr->ArrayLength >= arr->capacity)
            {
                arr->capacity = (arr->capacity * 3) / 2;
                arr->arr = safeReallocInt(arr->arr, arr->capacity);
            }

            if (readLong < INT_MIN || readLong > INT_MAX)
            {
                print("Integer overflow or underflow occurred\n");
                free(arr->arr);
                free(arr);
                exit(EXIT_FAILURE);
            }
            arr->arr[arr->ArrayLength++] = (int)readLong;
            ptr = endptr;
        }
    }

    return arr;
}
