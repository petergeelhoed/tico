#include "erf.h"
#include "mymath.h"

#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define BUFFER_SIZE 4096
#define ARR_SIZE 4096

int main(void)
{
    struct myarr* samples = readinput();

    struct stats result = fit_erf(samples->arrd, samples->ArrayLength);

    printf("σ=%lf μ=%lf N=%lu\n",
           result.stdev,
           result.mean,
           samples->ArrayLength);
}
