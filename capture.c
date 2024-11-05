#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <unistd.h>

#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"

#define ARR_BUFF 512

static int keepRunning = 1;
unsigned int columns = 80;
void sigint_handler(int signal)
{
    if (signal == SIGINT)
    {
        keepRunning = 0;
    }
    else if (signal == SIGWINCH)
    {
        struct winsize w;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
        columns = (unsigned int)w.ws_col;
        fprintf(stderr, "new width %d\n", columns);
    }
}

void set_signal_action(void)
{
    // Declare the sigaction structure
    struct sigaction act;

    // Set all of the structure's bits to 0 to avoid errors
    // relating to uninitialized variables...
    bzero(&act, sizeof(act));
    act.sa_handler = &sigint_handler;
    sigaction(SIGWINCH, &act, NULL);
    sigaction(SIGINT, &act, NULL);
}

int main(int argc, char* argv[])
{
    unsigned int rate = 48000;
    unsigned int bph = 21600;
    unsigned int evalue = 4; // width of gaussian window
    unsigned int zoom = 10;
    unsigned int time = 0;
    unsigned int everyline = 0;
    unsigned int len = 32;    //  syncwrite every len tics
    unsigned int cvalue = 7;  // cutoff for adding to correlation
    unsigned int verbose = 0; // print for this peak
    unsigned int fitN = 30;   // fit last 30 peaks, 10 seconds
    double SDthreshold = 3.;
    char* device = 0;
    FILE* rawfile = 0;
    FILE* fptotal = 0;
    FILE* fpDefPeak = 0;
    FILE* fpInput = 0;

    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int totalshift = 0;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

    int retVal = 0;
    int c;
    while ((c = getopt(argc, argv, "b:r:z:ht:s:e:c:d:w:p:f:D:v:I:l")) != -1)
    {
        switch (c)
        {
        case 'd':
            device = optarg;
            break;
        case 'f':
            retVal = checkUIntArg(c, &fitN, optarg);
            break;
        case 'c':
            retVal = checkUIntArg(c, &cvalue, optarg);
            cvalue = cvalue > 15 ? 15 : cvalue;
            break;
        case 'l':
            // number of characters to reserve for beaterror and rate
            everyline = 14;
            break;
        case 'v':
            retVal = checkUIntArg(c, &verbose, optarg);
            break;
        case 'e':
            retVal = checkUIntArg(c, &evalue, optarg);
            break;
        case 'I':
            retVal = checkFileArg(c, &fpInput, optarg, "r");
            break;
        case 'w':
            retVal = checkFileArg(c, &rawfile, optarg, "a");
            break;
        case 'D':
            retVal = checkFileArg(c, &fpDefPeak, optarg, "r");
            break;
        case 'p':
            retVal = checkFileArg(c, &fptotal, optarg, "w");
            break;
        case 's':
            SDthreshold = 0.0;
            SDthreshold = atof(optarg);
            if (SDthreshold == 0.0)
            {
                printf("invalid float argument for -s '%s'\n", optarg);
                return -1;
            }
            break;
        case 't':
            retVal = checkUIntArg(c, &time, optarg);
            break;
        case 'b':
            retVal = checkUIntArg(c, &bph, optarg);
            if (retVal == 0 && bph < 4800)
            {
                printf("refusing bph <4800 %d %d\n", bph, retVal);
                return -1;
            }
            break;
        case 'z':
            retVal = checkUIntArg(c, &zoom, optarg);
            break;
        case 'r':
            retVal = checkUIntArg(c, &rate, optarg);
            break;
        case 'h':
        default:
            fprintf(stderr,
                    "usage: capture \n"
                    "capture reads from the microphone and timegraphs your "
                    "mechanical watch\n"
                    "options:\n"
                    " -d default:1 capture device>\n"
                    " -z 10 zoom\n"
                    " -b 21600 bph of the watch\n"
                    " -r 48000 sampling rate in Hz\n"
                    " -t 30 measurment time in seconds\n"
                    " -s 3.0 cutoff standarddeviation\n"
                    " -w <file> write positions to file\n"
                    " -I <file> read from file instead of microphone\n"
                    " -p <file> write pulse to file\n"
                    " -D <file> read pulse from file\n"
                    " -c 7 threshold for local rate\n"
                    " -f 30 fit n points for local rate\n"
                    " -e 4 Gaussan convolution over input\n"
                    " -n 60 number of points to fit in local rate\n"
                    " -l print beaterror and rate on each line\n"
                    " -v <peak> write files for this peak \n");
            exit(0);
            break;
        }
        if (retVal != 0)
        {
            return retVal;
        }
    }

    // declarations
    unsigned int NN = rate * 7200 / bph;
    unsigned int mod = NN / zoom;
    unsigned int maxp = NN / 2;
    // should be even
    NN = (NN + NN % 2);
    unsigned int n = ARR_BUFF * 2;
    unsigned int tps = rate / NN;
    unsigned int maxtime = time ? time * tps : 30 * tps;
    fftw_complex* filterFFT = makeFilter(evalue, NN);
    int* maxpos = malloc(n * sizeof(int));
    int* maxvals = malloc(n * sizeof(int));
    int* derivative = malloc(NN * sizeof(int));
    int* reference = calloc(NN, sizeof(int));
    int* totaltick = calloc(NN, sizeof(int));

    device = (device == 0) ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = 0;
    if (fpInput == 0)
    {
        capture_handle = initAudio(format, device, rate);
    }
    char* buffer = malloc(NN * (unsigned int)snd_pcm_format_width(format) / 8);

    if (fpInput == 0 && capture_handle == 0)
    {
        fprintf(stderr, "No inputfile or soundcard");
        return -6;
    }
    if (buffer == 0 || totaltick == 0 || reference == 0 || maxvals == 0 ||
        maxpos == 0 || filterFFT == 0)
    {
        fprintf(stderr, "Could not allocate memory");
        return -5;
    }

    fprintf(stderr,
            "\033[2J\033[2;0H\nFound COLUMNS=%d, width = %.3fms  /  "
            "%.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));

    fillReference(fpDefPeak, reference, NN);

    unsigned int i = 0;
    unsigned int totalI = 0;
    while (keepRunning && !(totalI > maxtime && time))
    {
        if (i == n)
        {
            memcpy(maxpos, maxpos + ARR_BUFF, ARR_BUFF*sizeof(int));
            memcpy(maxvals, maxvals + ARR_BUFF, ARR_BUFF*sizeof(int));
            i -= ARR_BUFF;
        }

        int err = -32;
        while (err == -32)
        {
            int preshift = 0;

            if (i > 12)
            {
                preshift =
                    (((int)maxp + (int)NN / 2) % (int)(NN) - (int)(NN / 2));
                if (abs(preshift) > 1000)
                    preshift /= 8;
                else if (abs(preshift) > 100)
                    preshift /= 4;
                else if (abs(preshift) > 10)
                    preshift /= 2;
            }

            totalshift += preshift;

            err = readShiftedBuffer(
                derivative, capture_handle, NN, buffer, preshift, fpInput);
            if (err == -32)
            {
                fprintf(stderr, "Reinitializing capture_handle");
                if (rawfile)
                {
                    fprintf(rawfile, "# Reinitializing capture_handle");
                }
                snd_pcm_close(capture_handle);
                capture_handle = initAudio(format, device, rate);
                err = readBuffer(capture_handle, NN, buffer, derivative);
            }
            if (err == -33)
            {
                printf("Could not read integer from inputfile\n");
                return -7;
            }
        }
        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }
        if (totalI == 9)
        {
            // check after 8 ticktocks
            checkAndFlip(totaltick, reference, NN, verbose);
        }

        if (totalI == 12)
        {
            free(reference);
            reference = totaltick;
        }

        // NN/2 for no shift
        maxp = fftfit(derivative,
                      totaltick,
                      reference,
                      maxvals + i,
                      filterFFT,
                      NN,
                      totalI > 0 && totalI == verbose);

        maxpos[i] =
            totalshift + ((int)maxp + (int)NN / 2) % (int)(NN) - (int)(NN / 2);

        if (rawfile && i > 0 && i % len == 0)
        {
            syncappend(maxpos + i - len, len, rawfile);
        }

        fit10secs(&a, &b, &s, i, maxvals, maxpos, (int)cvalue, fitN);

        printheader(
            b * 86400 / NN, everyline, getBeatError(totaltick, NN, rate, 0));

        printspaces(maxpos[i],
                    maxvals[i],
                    (int)mod,
                    columns - everyline,
                    a,
                    (int)cvalue);
        i++;
        totalI++;
    }

    free(maxvals);
    free(buffer);
    free(derivative);
    fftw_free(filterFFT);

    if (rawfile)
    {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);

        struct tm* today = localtime(&tv.tv_sec);
        fprintf(rawfile,
                "# %04d-%02d-%02dT%02d:%02d:%02d.%ld %lu.%lu\n",
                today->tm_year + 1900,
                today->tm_mon + 1,
                today->tm_mday,
                today->tm_hour,
                today->tm_min,
                today->tm_sec,
                tv.tv_usec,
                tv.tv_sec,
                tv.tv_usec);

        writefile(rawfile, maxpos + i - i % len, i % len);
        fclose(rawfile);
    }
    free(maxpos);

    if (fptotal)
    {
        writefile(fptotal, totaltick, NN);
        fclose(fptotal);
    }

    if (capture_handle != 0)
    {
        snd_pcm_close(capture_handle);
    }

    free(totaltick);

    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    exit(0);
}
