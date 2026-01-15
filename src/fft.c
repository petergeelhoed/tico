#include <ctype.h>
#include <errno.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "myfft.h"
#include "mylib.h"
typedef struct
{
    unsigned int z;
    int pval;
    int lval;
    int fval;
} Config;

typedef struct
{
    double* data;
    unsigned int count;
} Signal;

double* safe_realloc(double* ptr, unsigned int new_size)
{
    double* next = realloc(ptr, new_size * sizeof(double));
    if (!next)
    {
        (void)fprintf(stderr, "Memory allocation failed\n");
        free(ptr);
        exit(-2);
    }
    return next;
}

Signal read_input(void)
{
    Signal signalStruct = {malloc(INIT_N * sizeof(double)), 0};
    unsigned int capacity = INIT_N;
    char line[BUFF_SIZE];
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

            if (signalStruct.count >= capacity)
            {
                capacity = (capacity * 3) / 2;
                signalStruct.data = safe_realloc(signalStruct.data, capacity);
            }

            signalStruct.data[signalStruct.count++] = readDouble;
            ptr = endptr;
        }
    }
    return signalStruct;
}

void run_fft(Signal sig, Config cfg)
{
    unsigned int arrayLength = sig.count * cfg.z;
    double par_a = 0.0;
    double par_b = 0.0;
    double s_err = 0.0;

    // Linear regression removal
    if (cfg.lval)
    {
        double* tmpx = malloc(sig.count * sizeof(double));
        for (unsigned int i = 0; i < sig.count; i++)
        {
            tmpx[i] = i;
        }
        linreg(tmpx, sig.data, sig.count, &par_a, &par_b, &s_err);
        (void)fprintf(
            stderr, "par_a=%lf par_b=%lf s=%lf\n", par_a, par_b, s_err);
        free(tmpx);
    }

    fftw_complex* data_in = fftw_alloc_complex(arrayLength);
    fftw_complex* data_out = fftw_alloc_complex(arrayLength);

    for (unsigned int i = 0; i < arrayLength; i++)
    {
        data_in[i][0] = (i < sig.count) ? (sig.data[i] - par_a - par_b * i) : 0;
        data_in[i][1] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(
        arrayLength, data_in, data_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Result Printing Logic
    if (cfg.pval)
    {
        for (unsigned int i = (cfg.lval ? 1 : 0); i < (sig.count + 2) / 2; i++)
        {
            double freq = cfg.fval ? (double)sig.count / (i * cfg.fval)
                                   : (double)i / sig.count;
            double mag = 2 * cfg.z *
                         sqrt(pow(data_out[i][0] / arrayLength, 2) +
                              pow(data_out[i][1] / arrayLength, 2));
            printf("%g %g\n", freq, mag);
        }
    }
    else
    {
        for (unsigned int i = 0; i < arrayLength; i++)
        {
            printf("%d %g %g\n",
                   i,
                   data_out[i][0] / arrayLength,
                   data_out[i][1] / arrayLength);
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(data_in);
    fftw_free(data_out);
}

int main(int argc, char** argv)
{
    Config cfg = {1, 0, 0, 0};
    int flag;
    int errno;
    char* endptr;

    while ((flag = getopt(argc, argv, "zplf:")) != -1)
    {
        switch (flag)
        {
        case 'z':
            cfg.z = 2;
            break;
        case 'p':
            cfg.pval = 1;
            break;
        case 'l':
            cfg.lval = 1;
            break;
        case 'f':
            errno = 0;
            endptr = NULL;
            long val = strtol(optarg, &endptr, DECIMAL);
            cfg.fval = (int)val;
            break;
        default:
            return -1;
        }
    }

    Signal sig = read_input();
    if (sig.count == 0)
    {
        free(sig.data);
        return -1;
    }

    run_fft(sig, cfg);

    free(sig.data);
    fftw_cleanup();
    return 0;
}
