#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"
#include "mymath.h"

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
    FILE* mfile = 0;
    FILE* fptotal = 0;
    FILE* fpDefPeak = 0;
    FILE* fpInput = 0;

    double b = 0.0;
    double a = 0.0;
    int totalshift = 0;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

    int retVal = 0;
    int c;
    while ((c = getopt(argc, argv, "b:r:z:ht:s:e:c:m:d:w:p:f:D:v:I:l")) != -1)
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
        case 'm':
            retVal = checkFileArg(c, &mfile, optarg, "w+");
            break;
        case 'w':
            retVal = checkFileArg(c, &rawfile, optarg, "w+");
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
                    " -m <file> write correlation to file\n"
                    " -I <file> read from file instead of microphone\n"
                    " -p <file> write pulse to file\n"
                    " -D <file> read pulse from file\n"
                    " -c 7 color crosscorrelation\n"
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
    unsigned int maxp = 0;
    // should be even
    NN = (NN + NN % 2);
    unsigned int n = ARR_BUFF * 2;
    unsigned int tps = rate / NN;
    unsigned int maxtime = time ? time * tps : 30 * tps;
    fftw_complex* filterFFT = makeFilter(evalue, NN);

    struct myarr maxpos = { malloc(n * sizeof(int)), 0 ,n};
    struct myarr maxvals = {0, calloc(n , sizeof(double)), n};

    struct myarr derivative = {malloc(NN * sizeof(int)), 0 ,NN}; 
    struct myarr reference = {malloc(NN * sizeof(int)), 0 ,NN};
    struct myarr totaltick = {malloc(NN * sizeof(int)), 0 ,NN};

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
    if (buffer == 0 || totaltick.arr == 0 || reference.arr == 0 || maxvals.arrd == 0 ||
        maxpos.arr == 0 || filterFFT == 0 || derivative.arr == 0)
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

    fillReference(fpDefPeak, &reference);

    unsigned int i = 0;
    unsigned int totalI = 0;
    while (keepRunning && !(totalI > maxtime && time))
    {
        if (i == n)
        {
            memcpy(maxpos.arr, maxpos.arr + ARR_BUFF, ARR_BUFF * sizeof(int));
            memcpy(maxvals.arrd, maxvals.arrd + ARR_BUFF, ARR_BUFF * sizeof(double));
            i -= ARR_BUFF;
        }

        int err = getData(maxp,
                          &totalshift,
                          rawfile,
                          fpInput,
                          capture_handle,
                          format,
                          device,
                          rate,
                          buffer,
                          derivative,
                          totalI);
        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }
        if (totalI == 9)
        {
            // check after 8 ticktocks
            checkAndFlip(&totaltick, &reference, verbose);
        }

        if (totalI == 12)
        {
            free(reference.arr);
            reference.arr = totaltick.arr;
        }

        // NN/2 for no shift
        maxp = fftfit(derivative,
                      totaltick.arr,
                      reference.arr,
                      maxvals.arrd + i,
                      filterFFT,
                      totalI > 0 && totalI == verbose);

        maxpos.arr[i] = totalshift + shiftHalf(maxp, NN);

        if (rawfile && i > 0 && i % len == 0)
        {
            syncappend(maxpos.arr + i - len, len, rawfile);
            syncappendDouble(maxvals.arrd + i - len, len, mfile);
        }

        fitNpeaks(&a, &b, i, &maxvals, &maxpos, fitN);

        printheader(
            b * 86400 / NN, everyline, getBeatError(&totaltick, rate, 0));

        printspaces(maxpos.arr[i],
                    (int)(maxvals.arrd[i]*16),
                    (int)mod,
                    columns - everyline,
                    a,
                    (int)cvalue);
        i++;
        totalI++;
    }

    free(buffer);
    free(derivative.arr);
    fftw_free(filterFFT);

    if (rawfile)
    {
        printTOD(rawfile);
        writefile(rawfile, maxpos.arr + i - i % len, i % len);
        if (mfile)
        {
            printTOD(mfile);
            writefileDouble(mfile, maxvals.arrd + i - i % len, i % len);
        }

        calculateTotalFromFile(totalI, rawfile, NN, SDthreshold);
        fclose(rawfile);
    }
    free(maxvals.arrd);
    free(maxpos.arr);

    if (fptotal)
    {
        writefile(fptotal, totaltick.arr, totaltick.NN);
        fclose(fptotal);
    }
    free(totaltick.arr);

    if (capture_handle != 0)
    {
        snd_pcm_close(capture_handle);
    }

    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    exit(0);
}
