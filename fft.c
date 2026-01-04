#include <errno.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "myfft.h"
#include "mylib.h"

#define INIT_N 4000
#define DECIMAL 10
#define BUFF_SIZE 256
int main(int argc, char** argv)
{
    unsigned int iN = INIT_N;
    unsigned int i = 0;

    int c;
    unsigned int z = 1;
    int pval = 0;
    int lval = 0;
    int fval = 0;
    char* endptr;

    while ((c = getopt(argc, argv, "zplf:")) != -1)
    {
        switch (c)
        {
        case 'l':
            lval = 1;
            break;
        case 'f':
            errno = 0; // Important: Clear errno before the call
            endptr = NULL;
            long val = strtol(optarg, &endptr, DECIMAL);
            fval = (int)val;
            break;
        case 'p':
            pval = 1;
            break;
        case 'z':
            z = 2;
            break;
        default:
            printf("usage\n fft -z zeropadding -p -l\n");
            printf("-p prints power spectrum i/length\n");
            printf("-l removes first order polynome\n");
            printf("-f frequency, will print 1/$1/freq\n");
            printf("period = length/$1/freq\n");
            printf("coefficient of sin = 2*sqrt($2*$2+$3*$3)\n");
            printf("l=63;s=2; seq 1 $((s*l)) | awk '{print sin($1/'$s')}'  | "
                   "fft -p | plot  'u (1/$1/'$s'):2 w l; set xrange [0:10]'\n");
            printf(" wav2raw < ../ruis.wav | awk '(NR>4096)'| fft -pl  |  plot "
                   " 'u ($1*44100):2 w l; set log'\n");
            return -1;
        }
    }

    double* tmpy = malloc(iN * sizeof(double));
    double* tmpx = malloc(iN * sizeof(double));
    double d;
    int ret = 0;
    char line[BUFF_SIZE];
    char* ptr;
    while (fgets(line, sizeof(line), stdin))
    {
        ptr = line;
        while (*ptr != '\0')
        {
            endptr = NULL;
            errno = 0;
            d = strtod(ptr, &endptr);
            if (ptr == endptr)
            {
                ptr++;
                continue; // No number found on this line
            }
            ptr = endptr;
            ret = 1;
            tmpy[i] = d;
            tmpx[i] = (double)i;

            i++;
            if (i == iN)
            {
                iN *= 3;
                iN /= 2;
                double* tmp2 = realloc(tmpy, iN * sizeof(double));
                if (tmp2)
                {
                    tmpy = tmp2;
                }
                else
                {
                    free(tmpx);
                    free(tmpy);
                    (void)fprintf(stderr, "Memory allocation failed");
                    return -2;
                }
                tmp2 = realloc(tmpx, iN * sizeof(double));
                if (tmp2)
                {
                    tmpx = tmp2;
                }
                else
                {
                    free(tmpx);
                    free(tmpy);
                    (void)fprintf(stderr, "Memory allocation failed");
                    return -2;
                }
            }
        }
    }
    if (ret == 0)
    {
        free(tmpx);
        free(tmpy);
        (void)fprintf(stderr, "Failed to parse double\n");
        return -1;
    }

    unsigned int N = i;
    double a = 0.0;
    double b = 0.0;
    double s = 0.0;
    if (lval)
    {
        linreg(tmpx, tmpy, N, &a, &b, &s);
        (void)fprintf(stderr, "a=%lf b=%lf s=%lf\n", a, b, s);
    }
    free(tmpx);

    fftw_complex* in = fftw_alloc_complex(N * z);
    for (i = 0; i < N; i++)
    {
        in[i][0] = tmpy[i] - a - b * i;
        in[i][1] = 0;
    }
    free(tmpy);

    for (i = N; i < N * z; i++)
    {
        in[i][0] = 0;
        in[i][1] = 0;
    }

    fftw_complex* out = fftw_alloc_complex(N * z);

    /* forward Fourier transform, save the result in 'out' */
    fftw_plan p =
        fftw_plan_dft_1d((int)(N * z), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    if (pval)
    {
        if (!lval)
        {
            printf("%f %g \n",
                   0.0,
                   sqrt(out[0][0] / (N * z) * out[0][0] / (N * z) +
                        out[0][1] / (z * N) * out[0][1] / (z * N)));
        }
        for (i = 1; i < (N + 2) / 2; i++)
        {
            printf("%g %g \n",
                   fval ? (double)(N) / (double)i / (double)fval
                        : (double)(i) / (double)N,
                   2 * z *
                       sqrt(out[i][0] / (N * z) * out[i][0] / (N * z) +
                            out[i][1] / (z * N) * out[i][1] / (z * N)));
        }
    }
    else
    {
        for (i = 0; i < N * z; i++)
        {
            printf("%d %g %g\n", i, out[i][0] / (N * z), out[i][1] / (z * N));
        }
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    return 0;
}
