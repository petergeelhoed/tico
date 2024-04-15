#include <math.h>
#include <fftw3.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mylib.h"

int main (int argc, char **argv)
{
    int N;
    int iN = 16000;
    fftw_complex *in, *out; /* double [2] */
    fftw_plan p;
    int i=0;

    int c;
    int z=1,pval=0,lval=0;

    while ((c = getopt (argc, argv, "zpl")) != -1)
    {
        switch (c)
        {
            case 'l':
                lval=1;
                break;
            case 'p':
                pval=1;
                break;
            case 'z':
                z=2;
                break;
            default:
                printf("usage\n fft -z zeropadding -p \n");
                printf("-p prints power spectrum length/i, coeff of sin\n");
                printf("period = Length/freq/$1\n");
                printf("coefficient of sin = 2*sqrt($2*$2+$3*$3)\n");
                return -1;
        }
    }

    float *tmp = malloc(iN*sizeof(float));
    float *tmpx = malloc(iN*sizeof(float));
    float d;
    while(scanf("%f ", &d) != EOF)
    {
        tmp[i] = d;
        tmpx[i] = i;

        i=i+1;
        if (i == iN) 
        {
            iN *= 2;
            float *tmp2 = realloc(tmp,iN*sizeof(float));
            if (tmp2)
            {
                tmp=tmp2;
            }
            tmp2 = realloc(tmpx,iN*sizeof(float));
            if (tmp2)
            {
                tmpx=tmp2;
            }
        }
    }


    N=i;
    double a=0,b=0,s=0;
    if (lval)
    {
        linregd(tmpx,tmp,N,&a,&b,&s);
        fprintf(stderr,"a=%lf b=%lf s=%lf\n",a,b,s);
    }



    in = fftw_alloc_complex(N*z);
    for (i = 0; i < N; i++) 
    {
        in[i][0] = tmp[i]-a-b*i;
        in[i][1] = 0;
    }

    for (i = N; i < N*z; i++) 
    {
        in[i][0] = 0;
        in[i][1] = 0;
    }

    out = fftw_alloc_complex(N*z);

    /* forward Fourier transform, save the result in 'out' */
    p = fftw_plan_dft_1d(N*z, in,out, FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(p);

    if (pval)
    {
        for (i = 1; i < (N+2)/2; i++) 
        {
            printf("%f %g \n",
                    //(float)i,
                    (float)(N)/i,
                    2*z*sqrt(
                        out[i][0]/(N*z)*out[i][0]/(N*z)+
                        out[i][1]/(z*N)*out[i][1]/(z*N))
                  ); 
        }
        if (!lval )
        {printf("%f %g \n",
                0.0,
                z*sqrt(
                    out[0][0]/(N*z)*out[0][0]/(N*z)+
                    out[0][1]/(z*N)*out[0][1]/(z*N))
              ); 
        }


    }else{
        for (i = 0; i < N*z; i++) 
        {
            printf("%3d ", i);
            printf("%g %g ", out[i][0]/(N*z), out[i][1]/(z*N) ); 
            printf("\n");
        }
    }
    fftw_destroy_plan(p);


    fftw_cleanup();
    return 0;
}
