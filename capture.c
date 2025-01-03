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
#include "mymath.h"
#include "mysound.h"
#include "mysync.h"

#define ARR_BUFF 512
#define DEFAULT_RATE 48000
#define DEFAULT_BPH 21600
#define DEFAULT_ZOOM 10
#define DEFAULT_EVALUE 4
#define DEFAULT_FITN 30
#define DEFAULT_LEN 32
#define DEFAULT_SDTHRESHOLD 3.0
#define DEFAULT_CVALUE 7
#define PRESHIFT_THRESHOLD 10
#define AUTOCOR_LIMIT 45
#define ERROR_SIGNAL -5

volatile int keepRunning = 1;
volatile unsigned int columns = 80;

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
    else
    {
        raise(signal);
    }
}

void set_signal_action(void)
{
    struct sigaction act;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    act.sa_handler = sigint_handler;

    // Set signal handlers
    if (sigaction(SIGWINCH, &act, NULL) == -1)
    {
        perror("sigaction");
        exit(ERROR_SIGNAL);
    }
    if (sigaction(SIGINT, &act, NULL) == -1)
    {
        perror("sigaction");
        exit(ERROR_SIGNAL);
    }
}

void parse_arguments(int argc,
                     char* argv[],
                     unsigned int* rate,
                     unsigned int* bph,
                     unsigned int* evalue,
                     unsigned int* zoom,
                     unsigned int* time,
                     unsigned int* everyline,
                     unsigned int* cvalue,
                     unsigned int* verbose,
                     unsigned int* fitN,
                     double* SDthreshold,
                     char** device,
                     FILE** rawfile,
                     FILE** mfile,
                     FILE** fptotal,
                     FILE** fpDefPeak,
                     FILE** fpInput)
{
    int c;
    while ((c = getopt(argc, argv, "b:r:z:ht:s:e:c:m:d:w:p:f:D:v:I:l")) != -1)
    {
        switch (c)
        {
        case 'd':
            *device = optarg;
            break;
        case 'f':
            if (checkUIntArg(c, fitN, optarg) != 0)
                exit(-1);
            break;
        case 'c':
            if (checkUIntArg(c, cvalue, optarg) != 0)
                exit(-1);
            *cvalue = *cvalue > 15 ? 15 : *cvalue;
            break;
        case 'l':
            *everyline = 14;
            break;
        case 'v':
            if (checkUIntArg(c, verbose, optarg) != 0)
                exit(-1);
            break;
        case 'e':
            if (checkUIntArg(c, evalue, optarg) != 0)
                exit(-1);
            break;
        case 'I':
            if (checkFileArg(c, fpInput, optarg, "r") != 0)
                exit(-1);
            break;
        case 'm':
            if (checkFileArg(c, mfile, optarg, "w+") != 0)
                exit(-1);
            break;
        case 'w':
            if (checkFileArg(c, rawfile, optarg, "w+") != 0)
                exit(-1);
            break;
        case 'D':
            if (checkFileArg(c, fpDefPeak, optarg, "r") != 0)
                exit(-1);
            break;
        case 'p':
            if (checkFileArg(c, fptotal, optarg, "w") != 0)
                exit(-1);
            break;
        case 's':
            *SDthreshold = atof(optarg);
            if (*SDthreshold == 0.0)
            {
                fprintf(stderr, "invalid float argument for -s '%s'\n", optarg);
                exit(-1);
            }
            break;
        case 't':
            if (checkUIntArg(c, time, optarg) != 0)
                exit(-1);
            break;
        case 'b':
            if (checkUIntArg(c, bph, optarg) != 0 || *bph < 4800)
            {
                fprintf(stderr, "refusing bph <4800 %d\n", *bph);
                exit(-1);
            }
            break;
        case 'z':
            if (checkUIntArg(c, zoom, optarg) != 0)
                exit(-1);
            break;
        case 'r':
            if (checkUIntArg(c, rate, optarg) != 0)
                exit(-1);
            break;
        case 'h':
        default:
            fprintf(stderr,
                    "usage: capture\n"
                    "capture reads from the microphone and timegraphs your "
                    "mechanical watch\n"
                    "options:\n"
                    " -d default:2 capture device>\n"
                    " -z 10 zoom\n"
                    " -b 21600 bph of the watch\n"
                    " -r 48000 sampling rate in Hz\n"
                    " -t 30 measurement time in seconds\n"
                    " -s 3.0 cutoff standard deviation\n"
                    " -w <file> write positions to file\n"
                    " -m <file> write correlation to file\n"
                    " -I <file> read from file instead of microphone\n"
                    " -p <file> write pulse to file\n"
                    " -D <file> read pulse from file\n"
                    " -c 7 cross-correlation limit\n"
                    " -f 30 fit n points for local rate\n"
                    " -e 4 Gaussian convolution over input\n"
                    " -n 60 number of points to fit in local rate\n"
                    " -l print beat error and rate on each line\n"
                    " -v <peak> write files for this peak\n");
            exit(0);
            break;
        }
    }
}

int main(int argc, char* argv[])
{
    unsigned int rate = DEFAULT_RATE;
    unsigned int bph = DEFAULT_BPH;
    unsigned int evalue = DEFAULT_EVALUE;
    unsigned int zoom = DEFAULT_ZOOM;
    unsigned int time = 0;
    unsigned int everyline = 0;
    unsigned int cvalue = DEFAULT_CVALUE;
    unsigned int verbose = 0;
    unsigned int len = DEFAULT_LEN;
    unsigned int fitN = DEFAULT_FITN;
    double SDthreshold = DEFAULT_SDTHRESHOLD;
    char* device = NULL;
    FILE* rawfile = NULL;
    FILE* mfile = NULL;
    FILE* fptotal = NULL;
    FILE* fpDefPeak = NULL;
    FILE* fpInput = NULL;

    parse_arguments(argc,
                    argv,
                    &rate,
                    &bph,
                    &evalue,
                    &zoom,
                    &time,
                    &everyline,
                    &cvalue,
                    &verbose,
                    &fitN,
                    &SDthreshold,
                    &device,
                    &rawfile,
                    &mfile,
                    &fptotal,
                    &fpDefPeak,
                    &fpInput);

    double a = 0.0;
    double b = 0.0;
    int totalshift = 0;
    struct winsize w;

    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

    unsigned int NN = rate * 7200 / bph;
    unsigned int maxp = 0;
    NN = (NN + NN % 2);
    unsigned int mod = NN / zoom;
    unsigned int n = ARR_BUFF * 2;
    unsigned int tps = rate / NN;
    unsigned int maxtime = time ? time * tps : 30 * tps;
    fftw_complex* filterFFT = makeFilter(evalue, NN);

    struct myarr maxpos = {malloc(n * sizeof(int)), 0, n};
    struct myarr maxvals = {0, calloc(n, sizeof(double)), n};
    struct myarr derivative = {malloc(NN * sizeof(int)), 0, NN};
    struct myarr tmpder = {malloc(NN * sizeof(int)), 0, derivative.NN};
    struct myarr reference = {malloc(NN * sizeof(int)), 0, NN};
    struct myarr totls[15];
    for (int t = 0; t < 15; t++)
    {
        totls[t].arr = malloc(NN * sizeof(int));
        totls[t].arrd = NULL;
        totls[t].NN = NN;
    }

    device = (device == NULL) ? "default:2" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = NULL;

    if (fpInput == NULL)
    {
        capture_handle = initAudio(format, device, rate);
    }
    char* buffer = malloc(NN * (unsigned int)snd_pcm_format_width(format) / 8);

    if (fpInput == NULL && capture_handle == NULL)
    {
        fprintf(stderr, "No inputfile or soundcard");
        return -6;
    }
    if (buffer == NULL || reference.arr == NULL || maxvals.arrd == NULL ||
        maxpos.arr == NULL || filterFFT == NULL || derivative.arr == NULL)
    {
        fprintf(stderr, "Could not allocate memory");
        return -5;
    }

    fprintf(stderr,
            "\033[2J\033[2;0H\nFound COLUMNS=%d, width = %.3fms / "
            "%.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));

    fillReference(fpDefPeak, &reference);

    struct sigaction sact;
    sigset_t new_set, old_set;

    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = sigint_handler;
    if (sigaction(SIGWINCH, &sact, NULL) != 0)
    {
        fprintf(stderr, "sigaction() error");
        return ERROR_SIGNAL;
    }
    sigemptyset(&new_set);
    sigaddset(&new_set, SIGWINCH);
    sigaddset(&new_set, SIGINT);

    unsigned int i = 0;
    unsigned int totalI = 0;
    while (keepRunning && !(totalI > maxtime && time))
    {
        if (i == n)
        {
            memcpy(maxpos.arr, maxpos.arr + ARR_BUFF, ARR_BUFF * sizeof(int));
            memcpy(maxvals.arrd,
                   maxvals.arrd + ARR_BUFF,
                   ARR_BUFF * sizeof(double));
            i -= ARR_BUFF;
        }

        // block signal to alsa
        if (sigprocmask(SIG_BLOCK, &new_set, &old_set) != 0)
        {
            fprintf(stderr, "sigprocmask() error");
            return ERROR_SIGNAL;
        }

        int err = getData(rawfile,
                          fpInput,
                          capture_handle,
                          format,
                          device,
                          rate,
                          buffer,
                          derivative);
        // restore signal handling
        if (sigprocmask(SIG_SETMASK, &old_set, NULL) != 0)
        {
            fprintf(stderr, "2nd sigprocmask() error");
            return ERROR_SIGNAL;
        }

        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }
        /*
                if (totalI == 9)
                {
                    checkAndFlip(&totaltick, &reference, verbose);
                }
                if (totalI == AUTOCOR_LIMIT)
                {
                    free(reference.arr);
                }
        */
        struct myarr* totaltick = &totls[totalI % 15];
        if (totalI > AUTOCOR_LIMIT)
            reference.arr = totaltick->arr;

        // totalshift modulo NN but that could also be negative
        // so modulate twice
        int totalm = (totalshift % (int)NN + (int)NN) % (int)NN;

        for (int j = 0; j < (int)NN; ++j)
        {
            // preshift the derivative
            tmpder.arr[j] = derivative.arr[(totalm + j) % (int)NN];
        }

        maxp = fftfit(tmpder,
                      totaltick->arr,
                      reference.arr,
                      maxvals.arrd + i,
                      filterFFT,
                      totalI > 0 && totalI == verbose);

        maxpos.arr[i] = totalshift + shiftHalf(maxp, NN);
        if (totalI > AUTOCOR_LIMIT && *(maxvals.arrd + i) > (double)cvalue / 16)
        {
            int preshift = shiftHalf(maxp, NN);

            if (abs(preshift) > PRESHIFT_THRESHOLD)
                preshift = (int)(3 * preshift / sqrt(abs(preshift)));

            totalshift += preshift;
        }

        if (rawfile && i > 0 && totalI % len == 0)
        {
            syncappend(maxpos.arr + i - len, len, rawfile);
            syncappendDouble(maxvals.arrd + i - len, len, mfile);
        }
        // should write all peaks, but this is just for diagnostics
        if (totalI % 15 == 0)
        {
            syncwrite(totaltick->arr, totaltick->NN, "livepeak");
        }

        fitNpeaks(&a, &b, i, &maxvals, &maxpos, fitN);

        printheader(b * 86400 / NN,
                    everyline,
                    getBeatError(totaltick, rate, 0),
                    (double)totalI / tps);
        printspaces(maxpos.arr[i],
                    (int)(maxvals.arrd[i] * 16),
                    (int)mod,
                    columns - everyline,
                    a,
                    (int)cvalue);

        i++;
        totalI++;
    }

    free(buffer);
    free(derivative.arr);
    free(tmpder.arr);
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

    for (int t = 0; t < 15; ++t)
    {
        struct myarr* totaltick = &totls[t];
        if (fptotal)
        {
            for (unsigned int j = 0; j < NN; ++j)
            {
                fprintf(fptotal, "%d %d %d\n", j, totaltick->arr[j], t);
            }
        }
        //       free(totaltick.arr);
    }
    if (fptotal)
    {
        fclose(fptotal);
    }

    if (capture_handle != NULL)
    {
        snd_pcm_close(capture_handle);
    }

    fprintf(stderr,
            "width = %.3fms / %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    return 0;
}
