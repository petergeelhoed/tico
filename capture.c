#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"

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
    unsigned int len = 30;    //  syncwrite every len tics
    unsigned int cvalue = 8;  // cutoff for adding to correlation
    unsigned int verbose = 0; // print for this peak
    unsigned int fitN = 30;   // fit last 30 peaks, 10 seconds
    double SDthreshold = 3.;
    char* device = 0;
    FILE* rawfile = 0;
    FILE* fptotal = 0;
    FILE* fpDefPeak = 0;
    FILE* fpInput = 0;

    int c;
    int retVal = 0;
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
            fprintf(
                stderr,
                "usage: capture \n"
                "capture reads from the microphone and timegraphs your watch\n"
                "options:\n"
                " -d <capture device> (default: 'default:1')\n"
                " -z <zoom> (default: 10)\n"
                " -b bph of the watch (default: 21600/h) \n"
                " -r sampling rate (default: 48000Hz)\n"
                " -t <measurment time> (default: 30s)\n"
                " -s cutoff standarddeviation (default: 3.0)\n"
                " -w <file> write positions to file\n"
                " -I <file> read from file instead of microphone\n"
                " -p <file> write pulse to file\n"
                " -D <file> read pulse from file\n"
                " -c 8 threshold for local rate\n"
                " -f 30 fit n points for local rate\n"
                " -e 4 Gaussan convolution over input\n"
                " -n 60 number of mpoints to fit in local rate\n"
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
    // should be even
    NN = (NN + NN % 2);
    unsigned int tps = rate / NN;
    unsigned int n = time ? time * tps : 30 * tps;
    unsigned int maxtime = n;
    int* maxpos = malloc(n * sizeof(int));
    int* maxvals = malloc(n * sizeof(int));
    unsigned int mod = NN / zoom;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    char spaces[1024];
    set_signal_action();

    device = device == 0 ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, rate);
    char* buffer = malloc(NN * (unsigned int)snd_pcm_format_width(format) / 8);

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    int* totaltick = malloc(NN * sizeof(int));
    for (unsigned int j = 0; j < NN; j++)
    {
        totaltick[j] = 0;
    }

    fprintf(stderr,
            "\033[2J\033[2;0H\nFound COLUMNS=%d, width = %.3fms  /  "
            "%.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));

    int* derivative = malloc(NN * sizeof(int));
    int* reference = malloc(NN * sizeof(int));

    fillReference(fpDefPeak, reference, NN);

    // read emptyparts
    readBuffer(capture_handle, NN, buffer, derivative);
    readBuffer(capture_handle, NN / 2, buffer, derivative);

    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int totalshift = 0;
    unsigned int maxp = 0;
    unsigned int i = 1;
    while (keepRunning && !(i > maxtime && time))
    {
        if (i == n)
        {
            n = n * 3 / 2;
            maxpos = realloc(maxpos, n * sizeof(int));
            maxvals = realloc(maxvals, n * sizeof(int));
        }
        int err = -32;
        while (err == -32)
        {
            err = readShiftedBuffer(derivative,
                                    capture_handle,
                                    NN,
                                    buffer,
                                    (int)maxp - (int)NN / 2,
                                    &totalshift,
                                    fpInput);
            if (err == -32)
            {
                snd_pcm_close(capture_handle);
                capture_handle = initAudio(format, device, rate);
                err = readBuffer(capture_handle, NN, buffer, derivative);
            }
        }

        if (i == 3 * tps)
        {
            checkAndFlip(totaltick, reference, NN, verbose);
        }

        if (i == 6 * tps)
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
                      i == verbose);

        if (rawfile && i > 0 && i % len == 0)
        {
            syncappend(maxpos + i - len, len, rawfile);
        }
        maxpos[i] = totalshift + ((int)maxp - (int)NN / 2);
        fit10secs(&a, &b, &s, i, maxvals, maxpos, (int)cvalue, fitN);
        printheader(
            b, NN, everyline, getBeatError(totaltick, NN, 0) * 1000. / rate);
        printspaces(maxpos[i],
                    maxvals[i],
                    spaces,
                    (int)mod,
                    columns - everyline,
                    a,
                    (int)cvalue);
        i++;
    }

    free(maxvals);
    free(buffer);
    fftw_free(filterFFT);
    snd_pcm_close(capture_handle);

    if (rawfile)
        syncappend(maxpos + i - i % len, i % len, rawfile);

    writefile(fptotal, totaltick, NN);
    if (fptotal)
    {
        fclose(fptotal);
    }

    calculateTotal(i, maxpos, NN, SDthreshold);
    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    free(maxpos);
    free(derivative);
    free(totaltick);
    exit(0);
}
