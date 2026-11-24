#include "myarr.h"
#include "printing.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define INIT_ARRSIZE 10
#define BUFFER_SIZE 512

void freemyarr(struct myarr* arr)
{
    free(arr->arr);
    free(arr->arrd);
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
    return ret;
}

/**
 * @brief Safely reallocates memory for a double array, exits on failure.
 *
 * @param ptr Pointer to the current memory block.
 * @param new_size New size for the array.
 * @return Pointer to the reallocated memory.
 */
static double* safeRealloc(double* ptr, unsigned int new_size)
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
struct myarr* readinput(void)
{
    unsigned int capacity = INIT_ARRSIZE;
    struct myarr* signalStruct = makemyarrd(capacity);
    char line[BUFFER_SIZE];
    char* ptr;
    char* endptr;
    size_t count = 0;

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

            if (count >= capacity)
            {
                capacity = (capacity * 3) / 2;
                signalStruct->arrd = safeRealloc(signalStruct->arrd, capacity);
            }

            signalStruct->arrd[count++] = readDouble;
            ptr = endptr;
        }
    }
    signalStruct->ArrayLength = count;

    return signalStruct;
}
