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

Signal read_input()
{
    Signal s = {malloc(INIT_N * sizeof(double)), 0};
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
            double d = strtod(ptr, &endptr);
            if (ptr == endptr)
            {
                ptr++;
                continue;
            }

            if (s.count >= capacity)
            {
                capacity = (capacity * 3) / 2;
                s.data = safe_realloc(s.data, capacity);
            }

            s.data[s.count++] = d;
            ptr = endptr;
        }
    }
    return s;
}

void run_fft(Signal sig, Config cfg)
{
    unsigned int Nz = sig.count * cfg.z;
    double a = 0.0;
    double b = 0.0;
    double s_err = 0.0;

    // Linear regression removal
    if (cfg.lval)
    {
        double* tmpx = malloc(sig.count * sizeof(double));
        for (unsigned int i = 0; i < sig.count; i++)
        {
            tmpx[i] = i;
        }
        linreg(tmpx, sig.data, sig.count, &a, &b, &s_err);
        (void)fprintf(stderr, "a=%lf b=%lf s=%lf\n", a, b, s_err);
        free(tmpx);
    }

    fftw_complex* in = fftw_alloc_complex(Nz);
    fftw_complex* out = fftw_alloc_complex(Nz);

    for (unsigned int i = 0; i < Nz; i++)
    {
        in[i][0] = (i < sig.count) ? (sig.data[i] - a - b * i) : 0;
        in[i][1] = 0;
    }

    fftw_plan p = fftw_plan_dft_1d(Nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Result Printing Logic
    if (cfg.pval)
    {
        for (unsigned int i = (cfg.lval ? 1 : 0); i < (sig.count + 2) / 2; i++)
        {
            double freq = cfg.fval ? (double)sig.count / (i * cfg.fval)
                                   : (double)i / sig.count;
            double mag = 2 * cfg.z *
                         sqrt(pow(out[i][0] / Nz, 2) + pow(out[i][1] / Nz, 2));
            printf("%g %g\n", freq, mag);
        }
    }
    else
    {
        for (unsigned int i = 0; i < Nz; i++)
        {
            printf("%d %g %g\n", i, out[i][0] / Nz, out[i][1] / Nz);
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

int main(int argc, char** argv)
{
    Config cfg = {1, 0, 0, 0};
    int c;
    int errno;
    char* endptr;

    while ((c = getopt(argc, argv, "zplf:")) != -1)
    {
        switch (c)
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
