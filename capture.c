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
    int cvalue = 5;
    int qvalue = 0;
    double threshold =3.;
    FILE* rawfile = 0;
    FILE* fptotal = fopen("total","w");
    if (fptotal == 0)
    {
        fprintf(stderr,"cannot open total\n");
        return -4;
    }

    while ((c = getopt (argc, argv, "b:r:z:ht:s:xwe:qc:d:")) != -1)
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
    int maxes[time * rate/buffer_frames];
    int maxvals[time * rate/buffer_frames];
    int mod = buffer_frames/mvalue;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);

    char *buffer = malloc(buffer_frames * snd_pcm_format_width(format) / 8);

    char out[16];
    int wdth;
    FILE *fp =popen("/usr/bin/tput cols" , "r");
    fgets(out,16,fp);
    fclose(fp);
    wdth=atoi(out);
    int columns = wdth - 10;
    char spaces[columns+1];

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

    int totaltick[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) totaltick[j] = 0;

    int totaltock[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) totaltock[j] = 0;

    int in[buffer_frames];
    int der[buffer_frames];
    unsigned char lsb;
    signed char msb;

    for (i = 0; i < time * rate/buffer_frames; ++i)
    {
        if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
        int max = 0;
        int maxpos = 0;
        for (int j = 0; j < buffer_frames*2; j+=2) {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2] = (msb << 8) | lsb ;

            int derivative = (j==0)?0:fabs(in[j/2]-in[j/2-1]);
            der[j/2] = derivative;
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
            int *reference = defaultpulse;
            int *total = (i%2==0||qvalue==0)?totaltick:totaltock;
            if (i>10*rate/buffer_frames)
            {
                reference = (i%2==0||qvalue==0)?totaltick:totaltock;
            }
            maxpos = fftfit(
                    der,
                    total,
                    reference,
                    &val,
                    filterFFT,
                    buffer_frames);
        }

        double b = 0.0;
        double a = 0.0;

        maxes[i] = maxpos;
        maxvals[i] = val;
        int n=0;
        int fitwindow = i>60?60:i;
        //int fitwindow = i/(1+qvalue);
        if (i >= fitwindow *(1+qvalue) )
        {
            int xarr[fitwindow];
            int yarr[fitwindow];
            double s = 0;
            for (int k = 0; k < fitwindow;k+=qvalue+1)
            {
                if (maxvals[i-k] > cvalue)
                {
                    yarr[n]=maxes[i-k];
                    xarr[n]=k;
                    n++;
                }
            }
            if (n > 1)
            {
                linreg(xarr,yarr, n, &a, &b, &s);
            }
        }

        int width = (maxpos%mod)*columns/mod;
        int widtha = (((int)a+mod)%mod)*columns/mod;
        fprintf(stderr,"%6.1fs/d",b*86400/buffer_frames);
        memset(spaces, ' ', columns);
        spaces[widtha] = '|';
        spaces[width] = '\0';
        fprintf(stderr,"%s%s%X\e[0m",spaces,i%2==0?"\e[31m": "\e[32m",val);
        memset(spaces, ' ', columns);
        if (widtha > width)
        {
            spaces[widtha-width-1] = '|';
            spaces[widtha-width-1+1] = '\0';
            fprintf(stderr,"%s",spaces);
        }
        fprintf(stderr,"\n");
    }

    for (int j = 0; j < buffer_frames; j++) fprintf(fptotal,"%d %d %d\n",totaltick[j],totaltock[j],defaultpulse[j]);

    free(buffer);

    snd_pcm_close (capture_handle);
    int n = time*rate/buffer_frames;

    double a = 0.0;
    double b = 0.0;
    double s = 0.0;
    int xarr[buffer_frames];
    for (int i = 0; i < n ; ++i) {xarr[i]=i;}

    linreg(xarr,maxes, n, &a, &b, &s);
    
    /*
    a /= buffer_frames*buffer_frames;
    b /= buffer_frames;
    s /= rate;
    */

    fprintf(stderr,"raw rate: %f s/d\n",-b*86400/buffer_frames);
    int m = 0;

    double e;

    for (int i = 0; i < n; ++i)
    {
        e = fabs(((double)maxes[i]-(a+xarr[i]*b))/s);
        if (e < threshold)
        {
            maxes[m] = maxes[i];
            xarr[m] = xarr[i];
            m++;
        }
    }
    linreg(xarr, maxes, m, &a, &b, &s);

    fprintf(stderr,"after %.1fσ removal: %.2f s/d\n",threshold,-b*86400/buffer_frames);
    fftw_free(filterFFT);
    exit (0);
}
