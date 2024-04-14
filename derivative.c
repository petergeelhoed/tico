#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void linreg(double* yarr, int NN, double* a, double* b, double* s , double *sa,double *sb)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (int j = -NN+1; j <= 0; ++j)
    {
        y  += yarr[j];
        xx += j*j;
        x  += j;
        xy += j*yarr[j];
        yy += yarr[j]*yarr[j];
    }

    *b = (NN*xy-x*y)/(NN*xx-x*x);
    *a = (y/NN)-((*b)*x/NN);
    *s = sqrt(fabs( 1./NN/(NN-2)*( NN*yy-y*y-(*b)*(*b)*(NN*xx-x*x))));
    
    *sb = sqrt(NN*(*s)/(NN*xx-x*x)) ;
    *sa = sqrt((*sb)*(*sb)/NN*xx);
}

int main (int argc, char *argv[])
{
    int bhp = 10;
    int pr = 1;
    int c;
    while ((c = getopt (argc, argv, "b:p:v")) != -1)
    {
        switch (c)
        {
            case 'v':
                fprintf (stderr,"x y a b s sa sb\n");
                break;
            case 'p':
                pr = atoi(optarg);
                break;
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
    double sa = 0.0;
    double sb = 0.0;
    int iN = 40;
    double *yarr = (double *)malloc(iN*sizeof(double));
    double d;
    int i =0; 
    int len = bhp;
    while(scanf("%lf", &d) != EOF)
    {
        yarr[i] = d;
        bhp = (len > i) ? i : len;
        if (i%pr==0)
        {
            if (i<3) 
            {
                a=yarr[0];
                b=d-a;
            }
            else
                linreg(yarr+i,bhp,&a, &b, &s,&sa,&sb);
            printf("%d %lf %lf %lf %lf %lf %lf\n",i,yarr[i],a,b,s,sa,sb);
        }
        i=i+1;
        if (i == iN) 
        {
            iN *= 2;
            yarr = (double *)realloc(yarr,iN*sizeof(double));
        }

    }

    free(yarr);
    exit (0);
}

