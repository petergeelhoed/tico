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
    int mvalue = 10;
    int time = 30;
    int c;
    char *device = 0;
    char defdev[200] = "default:1";
    int cvalue = 5;
    int qvalue = 0;
    double threshold =3.;
    FILE* rawfile = 0;
    FILE* fptotal = fopen("total","w");
    if (fptotal == 0)
    {
        fprintf(stderr,"cannot open file total\n");
        return -4;
    }

    while ((c = getopt (argc, argv, "b:r:z:ht:s:e:qc:d:w:")) != -1)
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
            case 'w':
                rawfile = fopen(optarg,"w");
                if (rawfile == 0)
                {
                    fprintf(stderr,"cannot open rawcapture\n");
                    return -4;
                }
                break;
            case 's':
                threshold = atof(optarg);
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
                        " -d <capture device> (default: default:1)\n"\
                        " -z <zoom> (default: 10)\n"\
                        " -b bph of the watch (default: 21600/h) \n"\
                        " -r sampling rate (default: 48000Hz)\n"\
                        " -t <measurment time> (default: 30s)\n"\
                        " -s cutoff standarddeviation (default: 3.0)\n"\
                        " -w <file> write positions to file"
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
    int NN = rate*3600/bph;
    int tps = rate/NN;
    device = device==0?defdev:device;

    fftw_complex *filterFFT = makeFilter(evalue, NN);
    int i;
    int n = time*tps; // total tics
    int maxpos[n];
    int maxvals[n];
    int mod = NN/mvalue;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);

    char *buffer = malloc(NN * snd_pcm_format_width(format) / 8);

    char out[16];
    FILE *fp =popen("/usr/bin/tput cols" , "r");
    fgets(out,16,fp);
    fclose(fp);
    int wdth=atoi(out);
    int columns = wdth - 10;
    char spaces[columns+1];

    fprintf(stderr,
            "Found COLUMNS=%d, width = %.3fms  /  %.1fμs/character\n",
            wdth - 1,
            mod*1000./rate,
            mod*1000000./rate/(wdth-1));

    int totaltick[NN];
    for (int j = 0; j < NN; j++) totaltick[j] = 0;

    int totaltock[NN];
    for (int j = 0; j < NN; j++) totaltock[j] = 0;


    // main loop
    int derivative[NN];
    int *reference = defaultpulse;
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    
    for (i = 0; i < n; ++i)
    {
        readBuffer(capture_handle, NN, buffer,derivative);
        if (i==10*tps) fprintf(stderr,"10 seconds, starting crosscor\n");

        if (i>10*tps)
        {
            reference = (i%2==0||qvalue==0)?totaltick:totaltock;
        }

        maxpos[i] = fftfit(
                derivative,
                (i%2==0||qvalue==0)?totaltick:totaltock,
                reference,
                maxvals+i,
                filterFFT,
                NN);

        fit10secs(&a,&b,&s,i,maxvals,maxpos,qvalue, cvalue);
        printspaces(maxpos[i],maxvals[i],spaces,mod,columns,a,b,NN,i);
    }

    free(buffer);
    fftw_free(filterFFT);
    snd_pcm_close (capture_handle);

    writefiles(fptotal, rawfile, totaltick, totaltick, defaultpulse, maxpos, n, NN);


    calculateTotal(n, maxpos, NN, threshold);
    exit (0);
}

