#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

#include "mylib.h"
#include "defaultpulse.h"

int main (int argc, char *argv[])
{

    unsigned int rate = 48000;
    int bph = 21600;
    int evalue = 4;
    int xvalue = 1;
    int mvalue = 10;
    int time = 30;
    int c;
    char *device = 0;
    char defdev[20] = "default:1";
    int cvalue = 8;
    int qvalue = 0;
    int verbose = 0;
    double threshold =3.;
    FILE* rawfile = 0;
    FILE* fptotal = fopen("total","w");
    if (fptotal == 0)
    {
        fprintf(stderr,"cannot open total\n");
        return -4;
    }

    while ((c = getopt (argc, argv, "b:r:z:ht:vs:xwe:qc:d:")) != -1)
    {
        switch (c)
        {
            case 'd':
                device = optarg;
                break;
            case 'c':
                cvalue = atoi(optarg);
                break;
            case 'q':
                qvalue = 1;
                break;
            case 'e':
                evalue = atoi(optarg);
                break;
            case 'x':
                xvalue = 0;
                break;
            case 'w':
                rawfile = fopen("rawcapture","w");
                if (rawfile == 0)
                {
                    fprintf(stderr,"cannot open rawcapture\n");
                    return -4;
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case 's':
                threshold = atoi(optarg);
                break;
            case 't':
                time = atoi(optarg);
                break;
            case 'b':
                bph = atoi(optarg);
                break;
            case 'z':
                mvalue = atoi(optarg);
                break;
            case 'r':
                rate = atoi(optarg);
                break;
            case 'h':
                fprintf (stderr,
                        "usage:\n"\
                        "capture " 
                        "options:\n"\
                        " -d <capture device> (default:1)\n"\
                        " -z <zoom> (default: 10)\n"\
                        " -b bph of the watch (default: 21600/h) \n"\
                        " -r sampling rate (default: 48000Hz)\n"\
                        " -t <measurment time> (default: 30s)\n"\
                        " -v verbose, print points to stdout\n"\
                        "    time, tick position modulo(3600/rate), deviation, σ\n"\
                        " -s cutoff standarddeviation (default: 3)\n"\
                        " -x do not use crosscorrelation instead use peak derivative\n"\
                        " -c 8 threshold for local rate\n"\
                        " -e 4 Gauss smooth\n"\
                        " -q split local tick/tock rate\n");
                exit(0);

            default:
                fprintf(stderr,"invalid option %c",c);
                exit(-2);
                break;
        }
    }
    int buffer_frames = rate*3600/bph;
    device = device==0?defdev:device;

    fftw_complex *filterFFT = makeFilter(evalue, buffer_frames);
    int i;
    int err;
    char *buffer;
    int maxes[time * rate/buffer_frames];
    int maxvals[time * rate/buffer_frames];
    int mod = buffer_frames/mvalue;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

    snd_pcm_t *capture_handle = initAudio(format, device, rate);


    buffer = malloc(buffer_frames * snd_pcm_format_width(format) / 8 );

    FILE *fp;
    char out[16];
    int wdth;
    fp =popen("/usr/bin/tput cols" , "r");
    fgets(out,16,fp);
    fclose(fp);
    wdth=atoi(out);

    fprintf(stderr,
            "Found COLUMNS=%d, width = %.3fms  /  %.1fμs/character\n",
            wdth - 1,
            mod*1000./rate,
            mod*1000000./rate/(wdth-1));

    for (i = 0; i < rate/buffer_frames; ++i)
    {
        if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
    }
    int total[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) total[j] = 0;

    for (i = 0; i < time * rate/buffer_frames; ++i)
    {
        if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
        unsigned char lsb;
        signed char msb;
        int in[buffer_frames];
        int der[buffer_frames];
        int max = 0;
        int maxpos = 0;
        for (int j = 0; j < buffer_frames*2; j+=2) {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2]=(msb << 8)| lsb ;

            int derivative = (j==0)?0:fabs(in[j/2]-in[j/2-1]);
            der[j/2] = derivative;
            //  printf("%d\n",derivative );
            if (derivative > max)
            {
                max = derivative;
                maxpos = j/2;
            }
        }

        if (i==10*rate/buffer_frames) fprintf(stderr,"10 seconds, starting crosscor\n");

        int val = 1;

        if (xvalue)
        {
            maxpos = fftfit(
                    der,
                    total,
                    i<10*rate/buffer_frames?defaultpulse:total,
                    &val,
                    filterFFT, buffer_frames);
        }

        double b=0.0;

        maxes[i] = maxpos;
        maxvals[i] = val;
        int n=0;
        if (i > 120)
        {
            double x = 0;
            double y = 0;
            double xx = 0;
            double xy = 0;
            double yy = 0;
            for (int k = 0; k < 60;k+=qvalue+1)
            {
                if (maxvals[i-k] > cvalue)
                {
                    n++;
                    y+=maxes[i-k];
                    xx+=k*k;
                    x+=k;
                    xy+=k*maxes[i-k];
                    yy+=maxes[i-k]*maxes[i-k];
                }
            }
            b = (n>1)?(n*xy-x*y)/(n*xx-x*x):0;


        }
        int columns = wdth - 1-9;
        int width = (maxpos%mod)*columns/mod;
        fprintf(stderr,"%6.1fs/d",b*86400/buffer_frames);
        char spaces[width+1];
        memset(spaces, ' ', width);
        spaces[width+1] = '\0';
        fprintf(stderr,"%s%s%X\e[0m\n",spaces,i%2==0?"\e[31m": "\e[32m",val);
    }

    for (int j = 0; j < buffer_frames; j++) fprintf(fptotal,"%d %d\n",total[j],defaultpulse[j]);

    free(buffer);

    snd_pcm_close (capture_handle);
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    int n = time*rate/buffer_frames;
    for (i = 0; i < n; ++i)
    {
        y+=maxes[i];
        xx+=i*i;
        x+=i;
        xy+=i*maxes[i];
        yy+=maxes[i]*maxes[i];
    }

    x=x*buffer_frames/rate;
    y=y/rate;
    yy=yy/rate/rate;
    xx=xx*buffer_frames/rate*buffer_frames/rate;
    xy=xy*buffer_frames/rate/rate;
    double a = (y*xx-x*xy)/(n*xx-x*x);
    double b=(n*xy-x*y)/(n*xx-x*x);
    double s = sqrt(( yy -2*a*y-2*b*xy+2*a*b*x+a*a*n+b*b*xx)/n);
    fprintf(stderr,"raw rate: %f s/d\n",-b*86400);
    x = 0.;
    y = 0.;
    xx = 0.;
    xy = 0.;
    yy = 0.;
    int m = 0;
    if (verbose)
    {
        for (i = 0; i < time * rate/buffer_frames; ++i)
        {
            printf("%f %f %f %f\n",
                    (double)i*buffer_frames/rate,(double)maxes[i]/rate,
                    ((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b)),s);
        }
    }

    double e;

    for (i = 0; i < n; ++i)
    {
        e = (((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
        if (e*e < threshold*threshold)
        {
            y+=maxes[i];
            xx+=i*i;
            x+=i;
            xy+=i*maxes[i];
            yy+=maxes[i]*maxes[i];
            m++;
        }
    }
    x=x*buffer_frames/rate;
    y=y/rate;
    yy=yy/rate/rate;
    xx=xx*buffer_frames/rate*buffer_frames/rate;
    xy=xy*buffer_frames/rate/rate;
    a = (y*xx-x*xy)/(m*xx-x*x);
    b=(m*xy-x*y)/(m*xx-x*x);
    s = sqrt(( yy -2*a*y-2*b*xy+2*a*b*x+a*a*m+b*b*xx)/m);
    fprintf(stderr,"after %.1fσ removal: %.2f s/d\n",threshold,-b*86400);
    /*  if (verbose)
        {
        for (i = 0; i < time * rate/buffer_frames; ++i)
        {
        e = (((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
        if (e*e < threshold*threshold)
        {
        printf("%f %f %f\n",
        (double)i*buffer_frames/rate,(double)maxes[i]/rate,
        ((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
        }
        }
        }
     */
    fftw_free(filterFFT);
    exit (0);
}
