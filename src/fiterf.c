#include "erf.h"

#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define BUFFER_SIZE 4096
#define ARR_SIZE 4096

static size_t readDoubles(double* samples)
{
    size_t size = 0;
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

            samples[size++] = readDouble;
            ptr = endptr;
        }
    }
    return size;
}
int main(void)
{
    double samples[ARR_SIZE];
    size_t size = readDoubles(samples);

    struct stats result = fit_erf(samples, size);

    printf("σ=%lf μ=%lf N=%lu\n", result.stdev, result.mean, size);
}
