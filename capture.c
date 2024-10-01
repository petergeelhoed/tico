#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "defaultpulse.h"
#include "myfft.h"
#include "mylib.h"
#include "mysound.h"

static int keepRunning = 1;
int columns = 80;
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
    int time = 0;
    int everyline = 0;
    int len = 30;     //  syncwrite every len tics
    int cvalue = 8;   // cutoff for adding to correlation
    int verbose = -1; // print for this peak
    int fitN = 30;    // fit last 30 peaks, 10 seconds
    double SDthreshold = 3.;
    char* device = 0;
    FILE* rawfile = 0;
    FILE* fptotal = 0;
    FILE* fpDefPeak = 0;
    FILE* fpInput = 0;

    int c;
    while ((c = getopt(argc, argv, "b:r:z:ht:s:e:c:d:w:p:f:D:v:I:l")) != -1)
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
            cvalue = cvalue > 15 ? 15 : cvalue;
            cvalue = cvalue < 0 ? 0 : cvalue;
            break;
        case 'l':
            everyline = 14;
            break;
        case 'v':
            verbose = atoi(optarg);
            break;
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
        case 'w':
            if (!access(optarg, F_OK))
            {
                fprintf(stderr, " existrawfile %s\n", optarg);
                if (remove(optarg))
                {
                    fprintf(stderr, "cannot delete rawfile\n");
                    return -4;
                }
            }
            rawfile = fopen(optarg, "a");
            if (rawfile == 0)
            {
                fprintf(stderr, "cannot open raw file for appending\n");
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
            SDthreshold = atof(optarg);
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
    }

    // declarations
    int NN = rate * 7200 / bph;
    // should be even
    NN = (NN + NN % 2);
    int tps = rate / NN;
    int n = time ? time * tps : 30 * tps;
    int maxtime = n;
    int* maxpos = malloc(n * sizeof(int));
    int* maxvals = malloc(n * sizeof(int));
    int mod = NN / zoom;

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    char spaces[1024];
    set_signal_action();

    device = device == 0 ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, rate);
    char* buffer = malloc(NN * snd_pcm_format_width(format) / 8);

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    int* totaltick = malloc(NN * sizeof(int));
    for (int j = 0; j < NN; j++)
        totaltick[j] = 0;

    fprintf(stderr,
            "\033[2J\033[2;0H\nFound COLUMNS=%d, width = %.3fms  /  "
            "%.1fμs/character\n",
            columns,
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));

    int* derivative = malloc(NN * sizeof(int));
    int* reference = malloc(NN * sizeof(int));
    int* defref = reference;

    // read default peak
    if (fpDefPeak != 0)
    {
        for (int j = 0; j < NN; j++)
        {

            if (fscanf(fpDefPeak, "%d", reference + j) != 1)
            {
                fprintf(stderr,
                        "not enough values in -D <default peak file>\n");

                exit(-5);
            }
        }
        fclose(fpDefPeak);
    }
    else if (NN == 16000)
    {
        reference = memcpy(reference, defaultpulsedouble, 16000 * sizeof(int));
    }
    else
    {
        for (int j = 0; j < NN; j++)
        {
            reference[j] = 0;
        }
        reference[NN / 4] = 1;
        reference[3 * NN / 4] = 1;
    }

    // read emptyparts
    readBuffer(capture_handle, NN, buffer, derivative);
    readBuffer(capture_handle, NN / 2, buffer, derivative);

    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int totalshift = 0;
    int maxp = 0;
    int i = 0;
    while (keepRunning && !(i > maxtime && time))
    {
        if (i == n)
        {
            n += 60 * tps;
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

        if (i == 3 * tps)
        {
            int* cross = malloc(NN * sizeof(int));
            crosscorint(NN, totaltick, reference, cross);
            int maxp = getmaxpos(cross, NN);
            if (verbose)
            {
                FILE* fp = fopen("flip", "w");
                for (int j = 0; j < NN; j++)
                {
                    fprintf(fp,
                            "%d %d %d %d\n",
                            j,
                            totaltick[j],
                            reference[j],
                            cross[j]);
                }
                fclose(fp);
            }
            if (maxp > NN / 4 && maxp < NN * 3 / 4)
            {
                fprintf(stderr, "FLIPPING peaks pos %d\n", maxp);

                int tmp = 0;
                for (int j = 0; j < NN / 2; j++)
                {
                    tmp = reference[j + NN / 2];
                    reference[j + NN / 2] = reference[j];
                    reference[j] = tmp;
                }
            }

            free(cross);
            readShiftedBuffer(derivative,
                              capture_handle,
                              NN,
                              buffer,
                              maxp,
                              &totalshift,
                              fpInput);
        }

        if (i == 6 * tps)
        {
            free(reference);
            defref = 0;
            reference = totaltick;
        }

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
        maxpos[i] = totalshift + maxp;
        fit10secs(&a, &b, &s, i, maxvals, maxpos, cvalue, fitN);
        printheader(
            b, NN, everyline, (getBeatError(totaltick, NN, 0)) * 1000. / rate);
        printspaces(
            maxpos[i], maxvals[i], spaces, mod, columns - everyline, a, cvalue);
        i++;
    }

    free(maxvals);
    free(buffer);
    fftw_free(filterFFT);
    snd_pcm_close(capture_handle);

    if (rawfile)
        syncappend(maxpos + i - i % len, i % len, rawfile);

    writefile(fptotal, totaltick, NN);
    fclose(fptotal);

    calculateTotal(i, maxpos, NN, SDthreshold);
    fprintf(stderr,
            "width = %.3fms  /  %.1fμs/character\n",
            mod * 1000. / rate,
            mod * 1000000. / rate / (columns - everyline));
    free(maxpos);
    free(derivative);
    free(totaltick);
    free(defref);
    exit(0);
}
