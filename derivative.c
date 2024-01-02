#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void linreg(double* xarr, double* yarr, int NN, double* a, double* b, double* s)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (int j = 0; j < NN; ++j)
    {
        y  += yarr[j];
        xx += xarr[j]*xarr[j];
        x  += xarr[j];
        xy += xarr[j]*yarr[j];
        yy += yarr[j]*yarr[j];
    }
    
    *a = (y*xx-x*xy)/(NN*xx-x*x);
    *b = (NN*xy-x*y)/(NN*xx-x*x);
    *s = sqrt((yy-2*(*a)*y-2*(*b)*xy+2*(*a)*(*b)*x+(*a)*(*a)*NN+(*b)*(*b)*xx)/NN);
}

int main (int argc, char *argv[])
{
    int bhp = 10;
    int c;
    while ((c = getopt (argc, argv, "b:")) != -1)
    {
        printf("%d %s %s\n", __LINE__,__func__,(__FILE__));
        switch (c)
        {
            case 'b':
                bhp = atoi(optarg);
                break;
            case 'h':
            default:
                fprintf (stderr,
                        "usage: derivative\n"\
                        " -b <b> option \n");
                exit(0);
                break;
        }
    }


    // declarations
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int iN = 40;
    double *yarr = (double *)malloc(iN*sizeof(double));
    double *xarr = (double *)malloc(iN*sizeof(double));
    double d;
    int i =0; 
    int len = bhp;
    while(scanf("%lf", &d) != EOF)
    {
        yarr[i] = d;
        xarr[i] = (double)i;
        bhp = (len > i) ? i : len;
        if (i>1)
            linreg(xarr+i-bhp,yarr+i-bhp,bhp,&a, &b, &s);
        printf("%lf %lf %lf %lf %lf\n",xarr[i],yarr[i],a,b,s);
        i=i+1;
        if (i == iN) 
        {
            iN *= 2;
            xarr = (double *)realloc(xarr,iN*sizeof(double));
            yarr = (double *)realloc(yarr,iN*sizeof(double));
        }

    }

    free(yarr);
    free(xarr);
    exit (0);
}

