#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mymath.h"
#include "mysignal.h"
#include "mysound.h"
#include "mysync.h"
#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#define KILO 1000.
#define MEGA 1000000.
#define ARR_BUFF 512
#define SECS_DAY 86400
#define HEXDEC 16
#define DEFAULT_RATE 48000
#define DEFAULT_BPH 21600
#define DEFAULT_ZOOM 10
#define DEFAULT_EVALUE 4
#define DEFAULT_FITN 30
#define DEFAULT_TIME 30
#define DEFAULT_TICKTOCK_WRITE 30
#define DEFAULT_TEETH 1
#define DEFAULT_SDTHRESHOLD 3.0
#define DEFAULT_CVALUE 7
#define PRESHIFT_THRESHOLD 100
#define PRESHIFT_THRESHOLD_ROOT 10
#define AUTOCOR_LIMIT 1
#define HALF 0.5
#define DEFAULT_COLUMNS 80
#define PCM_WIDTH 8
#define ERROR_ALLOCATE_MEM -5
#define ERROR_NO_SOURCE -6

volatile int keepRunning = 1;
volatile unsigned int columns = DEFAULT_COLUMNS;

// Struct to hold application resources for cleaner function signatures
typedef struct
{
    struct myarr* subpos;
    struct myarr* maxpos;
    struct myarr* maxvals;
    struct myarr* derivative;
    struct myarr* tmpder;
    struct myarr* reference;
    struct myarr** teethArray;
    fftw_complex* filterFFT;
    char* audioBuffer;
} AppResources;

void shift_buffer_data(unsigned int* ticktock,
                       struct myarr* subpos,
                       struct myarr* maxpos,
                       struct myarr* maxvals)
{
    memmove(subpos->arrd, subpos->arrd + ARR_BUFF, ARR_BUFF * sizeof(double));
    memmove(maxpos->arr, maxpos->arr + ARR_BUFF, ARR_BUFF * sizeof(int));
    memmove(maxvals->arrd, maxvals->arrd + ARR_BUFF, ARR_BUFF * sizeof(double));
    *ticktock -= ARR_BUFF;
}

int init_audio_source(CapConfig* cfg,
                      snd_pcm_t** handle,
                      unsigned int* actualRate)
{
    *actualRate = (unsigned int)(cfg->rate + HALF);
    if (cfg->fpInput == NULL)
    {
        printf("Casting inputrate %f to soundcard rate %d\n",
               cfg->rate,
               *actualRate);
        *handle = initAudio(SND_PCM_FORMAT_S16_LE, cfg->device, actualRate);
        printf("Actual rate %d, calculating with %f\n", *actualRate, cfg->rate);
    }
    return (*handle == NULL && cfg->fpInput == NULL) ? ERROR_NO_SOURCE : 0;
}

AppResources allocate_resources(unsigned int NN,
                                unsigned int ticktockBuffer,
                                unsigned int teeth,
                                int evalue)
{
    AppResources res;
    res.subpos = makemyarrd(ticktockBuffer);
    res.maxpos = makemyarr(ticktockBuffer);
    res.maxvals = makemyarrd(ticktockBuffer);
    res.derivative = makemyarr(NN);
    res.tmpder = makemyarr(NN);
    res.reference = makemyarr(NN);
    res.filterFFT = makeFilter(evalue, NN);
    res.audioBuffer = calloc(NN, 2); // 16-bit depth
    res.teethArray = malloc(sizeof(struct myarr*) * teeth);
    for (unsigned int t = 0; t < teeth; t++)
    {
        res.teethArray[t] = makemyarr(NN);
    }
    return res;
}

void cleanup_resources(AppResources* res, unsigned int teeth)
{
    free(res->audioBuffer);
    freemyarr(res->subpos);
    freemyarr(res->maxpos);
    freemyarr(res->maxvals);
    freemyarr(res->derivative);
    freemyarr(res->tmpder);
    freemyarr(res->reference);
    fftw_free(res->filterFFT);
    for (unsigned int t = 0; t < teeth; t++)
    {
        freemyarr(res->teethArray[t]);
    }
    free(res->teethArray);
}

void process_logging(CapConfig* cfg,
                     AppResources* res,
                     unsigned int tt,
                     unsigned int writeinterval)
{
    if (tt > 0 && tt % writeinterval == 0)
    {
        if (cfg->fpposition)
        {
            struct myarr* sync = makemyarrd(writeinterval);
            for (unsigned int k = 0; k < writeinterval; ++k)
            {
                sync->arrd[k] =
                    res->subpos->arrd[tt - writeinterval + k] +
                    (double)res->maxpos->arr[tt - writeinterval + k];
            }
            syncAppendMyarr(sync, cfg->fpposition);
            freemyarr(sync);
        }
        if (cfg->fpmaxcor)
        {
            struct myarr* sync = makemyarrd(writeinterval);
            memcpy(sync->arrd,
                   res->maxvals->arrd + tt - writeinterval,
                   writeinterval * sizeof(double));
            syncAppendMyarr(sync, cfg->fpmaxcor);
            freemyarr(sync);
        }
    }
}

// --- Main Program ---

int main(int argc, char* argv[])
{
    CapConfig cfg = {.rate = DEFAULT_RATE,
                     .bph = DEFAULT_BPH,
                     .evalue = DEFAULT_EVALUE,
                     .zoom = DEFAULT_ZOOM,
                     .fitN = DEFAULT_FITN,
                     .teeth = DEFAULT_TEETH,
                     .SDthreshold = DEFAULT_SDTHRESHOLD,
                     .cvalue = DEFAULT_CVALUE};
    parse_arguments(argc, argv, &cfg);

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

    snd_pcm_t* capture_handle = NULL;
    unsigned int actualRate;
    if (init_audio_source(&cfg, &capture_handle, &actualRate) != 0)
    {
        return ERROR_NO_SOURCE;
    }
    unsigned int NN = (actualRate * 2 * SECS_HOUR / cfg.bph);
    NN += (NN % 2);
    unsigned int mod = NN / cfg.zoom;
    const unsigned int maxtime =
        (unsigned int)cfg.rate * (cfg.time ? cfg.time : DEFAULT_TIME) / NN;

    AppResources res =
        allocate_resources(NN, ARR_BUFF * 2, cfg.teeth, cfg.evalue);
    fillReference(cfg.fpDefPeak, res.reference, cfg.teeth);

    sigset_t block;
    sigset_t non_block;
    setup_block_signals(&block);

    int totalshift = 0;
    int toothshift = 0;
    unsigned int ticktock = 0;
    unsigned int totalTickTock = 0;

    while (keepRunning && !(totalTickTock > maxtime && cfg.time))
    {
        if (ticktock == ARR_BUFF * 2)
        {
            shift_buffer_data(&ticktock, res.subpos, res.maxpos, res.maxvals);
        }
        block_signal(&block, &non_block);
        int err = getData(cfg.fpposition,
                          cfg.fpInput,
                          capture_handle,
                          SND_PCM_FORMAT_S16_LE,
                          cfg.device,
                          cfg.rate,
                          res.audioBuffer,
                          *res.derivative);
        unblock_signal(&non_block);
        if (err < 0)
        {
            break;
        }
        struct myarr* cumulativeTick =
            res.teethArray[totalTickTock % cfg.teeth];
        if (cfg.teeth > 1 && totalTickTock >= AUTOCOR_LIMIT * cfg.teeth)
        {
            toothshift = getshift(*res.teethArray[0],
                                  *res.teethArray[totalTickTock % cfg.teeth]);
        }
        for (unsigned int j = 0; j < NN; ++j)
        {
            res.tmpder->arr[j] =
                res.derivative->arr[modSigned(totalshift + j + toothshift, NN)];
        }
        int maxposition =
            shiftHalf(fftfit(*res.tmpder,
                             cumulativeTick->arr,
                             (totalTickTock < AUTOCOR_LIMIT * cfg.teeth)
                                 ? res.reference->arr
                                 : cumulativeTick->arr,
                             res.maxvals->arrd + ticktock,
                             res.filterFFT,
                             totalTickTock == cfg.verbose,
                             res.subpos->arrd + ticktock),
                      NN);

        res.maxpos->arr[ticktock] = totalshift + maxposition;
        if (totalTickTock > AUTOCOR_LIMIT &&
            res.maxvals->arrd[ticktock] > (double)cfg.cvalue / HEXDEC &&
            totalTickTock % cfg.teeth == 0)
        {
            if (abs(maxposition) > PRESHIFT_THRESHOLD)
            {
                maxposition = (int)(PRESHIFT_THRESHOLD_ROOT * maxposition /
                                    sqrt(abs(maxposition)));
            }
            totalshift += maxposition;
        }

        process_logging(&cfg, &res, ticktock, DEFAULT_TICKTOCK_WRITE);

        double a = 0.0;
        double b = 0.0;
        fitNpeaks(
            &a, &b, ticktock, res.maxvals, res.maxpos, res.subpos, cfg.fitN);
        printheader(b * SECS_DAY / NN,
                    cfg.everyline,
                    getBeatError(cumulativeTick, cfg.rate, 0),
                    (double)totalTickTock * NN / cfg.rate);
        printspaces(res.maxpos->arr[ticktock],
                    res.maxvals->arrd[ticktock] * HEXDEC,
                    mod,
                    columns - cfg.everyline,
                    a,
                    cfg.cvalue);

        ticktock++;
        totalTickTock++;
    }

    // Post-run cleanup and final file writes
    if (cfg.fpmaxcor)
    {
        printTOD(cfg.fpmaxcor);
        (void)fclose(cfg.fpmaxcor);
    }
    if (cfg.fpposition)
    {
        calculateTotalFromFile(
            totalTickTock, cfg.fpposition, NN, cfg.SDthreshold);
        (void)fclose(cfg.fpposition);
    }

    cleanup_resources(&res, cfg.teeth);
    if (capture_handle)
    {
        snd_pcm_close(capture_handle);
    }
    return 0;
}
