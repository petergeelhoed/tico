#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mymath.h"
#include "mysignal.h"
#include "mysound.h"
#include "mysync.h"

#define ARR_BUFF 512
#define DEFAULT_RATE 48000
#define DEFAULT_BPH 21600
#define DEFAULT_ZOOM 10
#define DEFAULT_EVALUE 4
#define DEFAULT_FITN 30
#define DEFAULT_TICKTOCK_WRITE 30
#define DEFAULT_TEETH 1
#define DEFAULT_SDTHRESHOLD 3.0
#define DEFAULT_CVALUE 7
#define PRESHIFT_THRESHOLD 10
#define AUTOCOR_LIMIT 1
#define ERROR_ALLOCATE_MEM -5
volatile int keepRunning = 1;
volatile unsigned int columns = 80;

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
                     unsigned int* teeth,
                     double* SDthreshold,
                     char** device,
                     FILE** fpposition,
                     FILE** fpmaxcor,
                     FILE** fptotal,
                     FILE** fpDefPeak,
                     FILE** fpInput)
{
    int c;
    while ((c = getopt(argc, argv, "b:r:z:ht:s:e:c:m:d:w:p:f:D:v:I:lj:")) != -1)
    {
        switch (c)
        {
        case 'd':
            *device = optarg;
            break;
        case 'j':
            if (checkUIntArg(c, teeth, optarg) != 0)
                exit(-1);
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
            if (checkFileArg(c, fpmaxcor, optarg, "w+") != 0)
                exit(-1);
            break;
        case 'w':
            if (checkFileArg(c, fpposition, optarg, "w+") != 0)
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
                    " -p <file> write pulse to file\n"
                    " -I <file> read from file instead of microphone\n"
                    " -D <file> read pulse from file\n"
                    " -c 7 cross-correlation limit\n"
                    " -f 30 fit n points for local rate\n"
                    " -e 4 Gaussian convolution over input\n"
                    " -n 60 number of points to fit in local rate\n"
                    " -l print beat error and rate on each line\n"
                    " -j 1 number of ratched wheel teeth\n"
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
    unsigned int writeinterval = DEFAULT_TICKTOCK_WRITE;
    unsigned int fitN = DEFAULT_FITN;
    unsigned int teeth = DEFAULT_TEETH;
    unsigned int ticktockBuffer = ARR_BUFF * 2;
    double SDthreshold = DEFAULT_SDTHRESHOLD;
    char* device = NULL;
    FILE* fpposition = NULL;
    FILE* fpmaxcor = NULL;
    FILE* fptotal = NULL;
    FILE* fpDefPeak = NULL;
    FILE* fpInput = NULL;
    double a = 0.0;
    double b = 0.0;
    int maxposition = 0;
    int totalshift = 0;
    struct winsize w;

    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

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
                    &teeth,
                    &SDthreshold,
                    &device,
                    &fpposition,
                    &fpmaxcor,
                    &fptotal,
                    &fpDefPeak,
                    &fpInput);

    unsigned int NN = rate * 7200 / bph;
    NN = (NN + NN % 2);
    unsigned int mod = NN / zoom;
    unsigned int tps = rate / NN;
    unsigned int maxtime = tps * (time ? time : 30);

    device = (device == NULL) ? "default:2" : device;

    // initialize sound source
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = NULL;

    if (fpInput == NULL)
    {
        capture_handle = initAudio(format, device, rate);
    }

    if (fpInput == NULL && capture_handle == NULL)
    {
        fprintf(stderr, "No inputfile or soundcard");
        return -6;
    }

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    struct myarr maxpos = {
        calloc(ticktockBuffer, sizeof(int)), 0, ticktockBuffer};
    struct myarr maxvals = {
        0, calloc(ticktockBuffer, sizeof(double)), ticktockBuffer};
    struct myarr derivative = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr tmpder = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr reference = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr teethArray[teeth];
    for (unsigned int t = 0; t < teeth; t++)
    {
        teethArray[t].arr = calloc(NN, sizeof(int));
        if (teethArray[t].arr == NULL)
        {
            fprintf(stderr, "Could not allocate memory");
            return ERROR_ALLOCATE_MEM;
        }

        teethArray[t].arrd = NULL;
        teethArray[t].NN = NN;
    }

    char* buffer = calloc(NN, (unsigned int)snd_pcm_format_width(format) / 8);
    if (buffer == NULL || reference.arr == NULL || maxvals.arrd == NULL ||
        maxpos.arr == NULL || filterFFT == NULL || derivative.arr == NULL ||
        tmpder.arr == NULL)
    {
        fprintf(stderr, "Could not allocate memory");
        return ERROR_ALLOCATE_MEM;
    }

    fprintf(stderr,
            "\033[2J\033[2;0H\nFound COLUMNS=%d, width = %.3fms / "
            "%.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));

    fillReference(fpDefPeak, &reference, teeth);

    sigset_t block;
    sigset_t non_block;
    setup_block_signals(&block);
    int toothshift = 0;
    unsigned int ticktock = 0;
    unsigned int totalTickTock = 0;
    while (keepRunning && !(totalTickTock > maxtime && time))
    {
        if (ticktock == ticktockBuffer)
        {
            // shift data back, has been written already
            memcpy(maxpos.arr, maxpos.arr + ARR_BUFF, ARR_BUFF * sizeof(int));
            memcpy(maxvals.arrd,
                   maxvals.arrd + ARR_BUFF,
                   ARR_BUFF * sizeof(double));
            ticktock -= ARR_BUFF;
        }

        block_signal(&block, &non_block);
        int err = getData(fpposition,
                          fpInput,
                          capture_handle,
                          format,
                          device,
                          rate,
                          buffer,
                          derivative);
        unblock_signal(&non_block);

        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }

        struct myarr* cumulativeTick = &teethArray[totalTickTock % teeth];
        if (totalTickTock >= AUTOCOR_LIMIT * teeth)
        {
            // make sure this is only done after the j teethArray are filled at
            // least once
            if (teeth > 1)
            {
                toothshift =
                    getshift(teethArray[0], teethArray[totalTickTock % teeth]);
            }

            if (totalTickTock == AUTOCOR_LIMIT * teeth)
            {
                free(reference.arr);
            }
            // use the appropriate tick as a reference
            reference.arr = cumulativeTick->arr;
        }

        for (unsigned int j = 0; j < NN; ++j)
        {
            // preshift the derivative
            tmpder.arr[j] =
                derivative.arr[modSigned(totalshift + j + toothshift, NN)];
        }

        maxposition =
            shiftHalf(fftfit(tmpder,
                             cumulativeTick->arr,
                             reference.arr,
                             maxvals.arrd + ticktock,
                             filterFFT,
                             totalTickTock > 0 && totalTickTock == verbose),
                      NN);

        maxpos.arr[ticktock] = totalshift + maxposition;

        if (totalTickTock > AUTOCOR_LIMIT &&
            *(maxvals.arrd + ticktock) > (double)cvalue / 16)
        {
            if (abs(maxposition) > PRESHIFT_THRESHOLD)
            {
                // shift with 3√ outside threshold
                maxposition = (int)(3 * maxposition / sqrt(abs(maxposition)));
            }

            totalshift += maxposition;
        }

        if (ticktock > 0 && totalTickTock % writeinterval == 0)

        {
            if (fpposition)
            {
                syncappend(maxpos.arr + ticktock - writeinterval,
                           writeinterval,
                           fpposition);
            }
            if (fpmaxcor)
            {
                syncappendDouble(maxvals.arrd + ticktock - writeinterval,
                                 writeinterval,
                                 fpmaxcor);
            }
        }

        fitNpeaks(&a, &b, ticktock, &maxvals, &maxpos, fitN);

        printheader(b * 86400 / NN,
                    everyline,
                    getBeatError(cumulativeTick, rate, 0),
                    (double)totalTickTock / tps);
        printspaces(maxpos.arr[ticktock],
                    maxvals.arrd[ticktock] * 16,
                    mod,
                    columns - everyline,
                    a,
                    cvalue);

        ticktock++;
        totalTickTock++;
    }

    if (teeth > 1)
    {
        printf("peak   shift beaterr\n");
        for (unsigned int k = 0; k < teeth; ++k)
        {
            double beaterr = getBeatError(&teethArray[k], rate, 0);
            printf("%6d%6d%6.2f\n",
                   k,
                   getshift(teethArray[0], teethArray[k]),
                   beaterr);
        }
    }

    free(buffer);
    free(derivative.arr);
    free(tmpder.arr);
    fftw_free(filterFFT);

    thread_lock();
    if (fpmaxcor)
    {
        printTOD(fpmaxcor);
        writefileDouble(fpmaxcor,
                        maxvals.arrd + ticktock - ticktock % writeinterval,
                        ticktock % writeinterval);
    }
    if (fpposition)
    {
        printTOD(fpposition);
        writefile(fpposition,
                  maxpos.arr + ticktock - ticktock % writeinterval,
                  ticktock % writeinterval);
        //   syncappend(maxpos.arr + ticktock - ticktock % writeinterval,
        //   ticktock % writeinterval, fpposition);
        calculateTotalFromFile(totalTickTock, fpposition, NN, SDthreshold);
        fclose(fpposition);
    }
    thread_unlock();

    free(maxvals.arrd);
    free(maxpos.arr);

    for (unsigned int t = 0; t < teeth; ++t)
    {
        struct myarr* cumulativeTick = &teethArray[t];
        if (fptotal)
        {
            int toothshift = getshift(teethArray[0], *cumulativeTick);
            for (unsigned int j = 0; j < NN; ++j)
            {
                fprintf(fptotal,
                        "%d %d %u %d\n",
                        shiftHalf(j + toothshift, NN),
                        cumulativeTick->arr[j],
                        t,
                        toothshift);
            }
        }
    }
    for (unsigned int t = 0; t < teeth; ++t)
    {
        free(teethArray[t].arr);
    }
    if (fpInput)
    {
        fclose(fpInput);
    }
    if (fptotal)
    {
        fclose(fptotal);
    }

    if (capture_handle != NULL)
    {
        snd_pcm_close(capture_handle);
        snd_pcm_hw_free(capture_handle);
    }

    fprintf(stderr,
            "width = %.3fms / %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    return 0;
}
