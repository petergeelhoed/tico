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
    int kvalue = 1;
    int zoom = 10;
    int time = 30;
    int c;
    int cvalue = 5;
    int vvalue = -1;
    int fitN = 60;
    int qvalue = 0;
    char *device = 0;
    double threshold =3.;
    FILE* rawfile = 0;
    FILE* fptotal = 0;
    FILE* fpDefPeak = 0;

    while ((c = getopt (argc, argv, "b:r:z:ht:s:e:qc:d:w:p:f:kD:v:")) != -1)
    {
        switch (c)
        {
            case 'd':
                device = optarg;
                break;
            case 'f':
                fitN = atoi(optarg);
                break;
            case 'c':
                cvalue = atoi(optarg);
                break;
            case 'v':
                vvalue = 1;
                break;
            case 'q':
                qvalue = 1;
                break;
            case 'k':
                kvalue = 0;
                break;
            case 'e':
                evalue = atoi(optarg);
                break;
            case 'w':
                rawfile = fopen(optarg, "w");
                if (rawfile == 0)
                {
                    fprintf(stderr, "cannot open rawcapture\n");
                    return -4;
                }
                break;
            case 'D':
                fpDefPeak = fopen(optarg, "r");
                if (fpDefPeak == 0)
                {
                    fprintf(stderr, "cannot open file -D <file>\n");
                    return -4;
                }
                break;
            case 'p':
                fptotal = fopen(optarg, "w");
                if (fptotal == 0)
                {
                    fprintf(stderr, "cannot open file -p <file>\n");
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
                zoom = atoi(optarg);
                break;
            case 'r':
                rate = atoi(optarg);
                break;
            case 'h':
            default:
                fprintf (stderr,
                        "usage: capture \n"\
                        "capture reads from the microphone and timegraphs your watch\n" 
                        "options:\n"\
                        " -d <capture device> (default: 'default:1')\n"\
                        " -z <zoom> (default: 10)\n"\
                        " -b bph of the watch (default: 21600/h) \n"\
                        " -r sampling rate (default: 48000Hz)\n"\
                        " -t <measurment time> (default: 30s)\n"\
                        " -s cutoff standarddeviation (default: 3.0)\n"\
                        " -w <file> write positions to file\n"
                        " -p <file> write pulse to file\n"
                        " -D <file> read pulse from file\n"
                        " -c 8 threshold for local rate\n"\
                        " -e 4 Gauss smooth\n"\
                        " -k do not correlate tick and tock together\n"\
                        " -n 60 number of mpoints to fit in local rate\n"\
                        " -v <peak> write files for this peak \n"\
                        " -q split local tick/tock rate\n");
                exit(0);
                break;
        }
    }

    if (qvalue && kvalue)
    {
        fprintf (stderr,"cannot ue -q and -k together\n");
        exit (-3);

    }

    // declarations
    int NN = rate*3600/bph*(kvalue+1);
    // should be even
    NN = (NN+NN%2);

    int tps = rate/NN;
    int n = time*tps; 
    int maxpos[n];
    int maxvals[n];
    int mod = NN/zoom;
    
    FILE *fp =popen("/usr/bin/tput cols", "r");
    char out[16];
    fgets(out, 16, fp);
    fclose(fp);

    int wdth=atoi(out);
    int columns = wdth - 17;
    char spaces[columns+1];

    device = device==0?"default:1":device;

    fftw_complex *filterFFT = makeFilter(evalue, NN);
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);
    char *buffer = malloc(NN * snd_pcm_format_width(format) / 8);
    
    int* totaltick = malloc(NN*sizeof(int));
    for (int j = 0; j < NN; j++) totaltick[j] = 0;

    int* totaltock = malloc(NN*sizeof(int));
    for (int j = 0; j < NN; j++) totaltock[j] = 0;



    fprintf(stderr,
            "Found COLUMNS=%d, width = %.3fms  /  %.1fμs/character\n",
            wdth - 1,
            mod*1000./rate,
            mod*1000000./rate/(wdth-1));

    // main loop
    int *derivative = malloc(NN*sizeof(int));
    int *reference;
    int *referenceBase;

    if (fpDefPeak != 0)
    {
        reference = malloc(NN*sizeof(int));
        for (int j = 0; j < NN; j++) 
        {

            if (fscanf(fpDefPeak,"%d",reference+j) != 1)
            {
                fprintf(stderr, "not enough values in -D <default peak file>\n");
                
                exit(-5);
            }

        }

        referenceBase = reference;
    }
    else if (NN==8000)
    {
        reference = defaultpulse;
        referenceBase = reference;
    }
    else if (NN==16000 && kvalue == 1)
    {
        reference = defaultpulsedouble;
        referenceBase = reference;

 /*       reference = malloc(NN*sizeof(int));
        memcpy(reference,defaultpulse,8000*sizeof(int));
        memcpy(reference+8000,defaultpulse,8000*sizeof(int));
        referenceBase = malloc(NN*sizeof(int));
        memcpy(referenceBase,reference,16000*sizeof(int));
        */
    }
    else
    {
        reference = malloc(NN*sizeof(int));

        for (int j=0;j<NN;j++)
        {
            reference[j] = 0;
        }

        int min = (NN-8000)/2;
        if (min > 0)
        {
            memcpy(reference+min,defaultpulse,8000*sizeof(int));
        }
        else
        {
            memcpy(reference,defaultpulse-min,(8000+2*(min))*sizeof(int));
        }
    }

    // read emptyparts 
    readBuffer(capture_handle, NN, buffer, derivative);
    readBuffer(capture_handle, NN/2, buffer, derivative);

    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int i = 0;
    int totalshift = 0;
    int bound = 32;
    int upperBound = +NN/bound;
    int lowerBound = -NN/bound;
    int shift = NN/bound/10;
    int maxp = 0;
    for (; i < n; ++i)
    {

        readShiftedBuffer(derivative, capture_handle, NN, buffer, maxp, shift, &totalshift, lowerBound, upperBound);

        if (i==3*tps)
        {
            int* cross = malloc(NN*sizeof(int));
            crosscorint(NN, totaltick, reference,cross);
            int maxp = 0;
            int maxval = -NN*NN;
            for (int j = 0; j < NN; j++)
            {
                if (maxval < cross[j])
                {
                    maxp= j;
                    maxval = cross[j];
                }
            }
            if (maxp > NN/4 && maxp < NN*3/4)
            {
                fprintf(stderr,"FLIPPING peaks pos %d val %d\n",maxp,maxval);

                int tmp =0;
                for (int j = 0; j < NN/2; j++)
                {
                    tmp = reference[j+NN/2];
                    reference[j+NN/2] = reference[j];
                    reference[j] = tmp;
                }

            }

            free(cross);

        }
        

        if (i>6*tps)
        {
            reference = (i%2==0||qvalue==0)?totaltick:totaltock;
        }

        maxp = fftfit(
                derivative,
                (i%2==0||qvalue==0)?totaltick:totaltock,
                reference, maxvals+i, filterFFT, NN, kvalue,
                i==vvalue);

        maxpos[i] = totalshift + maxp;


        fit10secs(&a, &b, &s, i, maxvals, maxpos, qvalue, cvalue, fitN);
        printspaces(maxpos[i], maxvals[i], spaces, mod, columns, a, b, NN, kvalue?0:i, (float)(getBeatError(totaltick, NN,i==vvalue))/rate*1000);
    }

    free(buffer);
    fftw_free(filterFFT);
    snd_pcm_close (capture_handle);

    writefiles(fptotal, rawfile, totaltick, totaltock, referenceBase, maxpos, n, NN);

    calculateTotal(n, maxpos, NN, threshold);
    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod*1000./rate,
            mod*1000000./rate/(wdth-1));
    free(derivative);
    free(totaltick);
    free(totaltock);
    free(reference);
    exit (0);
}

