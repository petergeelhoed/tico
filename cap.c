#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>

#include "myfft.h"
#include "mylib.h"
#include "mysound.h"

static int keepRunning = 1;
static int columns = 80;
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
        columns = w.ws_col;
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
    int bph = 21600;
    int evalue = 4; // width of gaussian window
    int zoom = 10;
    int time = 5;
    int cvalue = 8; // cutoff for adding to correlation
    int fitN = 30;  // fit last 30 peaks, 10 seconds
    double SDthreshold = 3.;
    char device[] = "default:1";
    FILE* fpInput = 0;

    int c;
    while ((c = getopt(argc, argv, "e:z:I:h")) != -1)
    {
        switch (c)
        {
        case 'e':
            evalue = atoi(optarg);
            break;
        case 'I':
            fpInput = fopen(optarg, "r");
            if (fpInput == 0)
            {
                fprintf(stderr, "cannot open input file\n");
                return -4;
            }
            break;
        case 'z':
            zoom = atoi(optarg);
            break;
        case 'h':
        default:
            fprintf(
                stderr,
                "usage: capture \n"
                "capture reads from the microphone and timegraphs your watch\n"
                "options:\n"
                " -z <zoom> (default: 10)\n"
                " -I <file> read from file instead of microphone\n"
                " -e 4 Gaussan convolution over input\n");
            exit(0);
            break;
        }
    }

    // declarations
    int NN = rate * 7200 / bph;
    // should be even
    NN = (NN + NN % 2);
    int tps = rate / NN;
    int n = time * tps;
    int* maxpos = malloc(n * sizeof(int));
    int* maxvals = malloc(n * sizeof(int));
    int mod = NN / zoom;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    char spaces[1024];
    set_signal_action();

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, rate);
    char* buffer = malloc(NN * snd_pcm_format_width(format) / 8);

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    int* totaltick = malloc(NN * sizeof(int));
    for (int j = 0; j < NN; j++)
    {
        totaltick[j] = 0;
    }
    totaltick[NN / 4] = 1;
    totaltick[3 * NN / 4] = 1;

    fprintf(stderr,
            "Found COLUMNS=%d, width = %.3fms  /  %.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - 14));

    int* derivative = malloc(NN * sizeof(int));

    // read emptyparts
    readBuffer(capture_handle, NN, buffer, derivative);
    readBuffer(capture_handle, NN / 2, buffer, derivative);

    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int totalshift = 0;
    int maxp = 0;
    int i = 0;
    while (keepRunning)
    {
        if (i == n)
        {
            n *= 1.5;
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
                                    maxp,
                                    &totalshift,
                                    fpInput);
            if (err == -32)
            {
                snd_pcm_close(capture_handle);
                capture_handle = initAudio(format, device, rate);
                err = readBuffer(capture_handle, NN, buffer, derivative);
            }
        }

        maxp = fftfit(
            derivative, totaltick, totaltick, maxvals + i, filterFFT, NN, 0);

        maxpos[i] = totalshift + maxp;
        fit10secs(&a, &b, &s, i, maxvals, maxpos, cvalue, fitN);
        printheader(b, NN, 0, (getBeatError(totaltick, NN, 0)) * 1000. / rate);
        printspaces(
            maxpos[i], maxvals[i], spaces, mod, columns - 14, a, cvalue);
        i++;
    }

    free(maxvals);
    free(buffer);
    fftw_free(filterFFT);
    snd_pcm_close(capture_handle);

    calculateTotal(i, maxpos, NN, SDthreshold);
    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - 14));
    free(maxpos);
    free(derivative);
    free(totaltick);
    exit(0);
}
