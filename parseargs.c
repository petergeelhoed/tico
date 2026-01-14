#include "parseargs.h"
#include "mylib.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void print_usage(void)
{
    (void)fprintf(stderr,
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
}
int checkFloatArg(int name, double* value, char* optarg)
{

    char* endp = NULL;
    *value = strtod(optarg, &endp);
    if (*value == 0 || optarg == endp)
    {
        printf("invalid float argument for -%c: '%s'\n", (char)name, optarg);
        return -1;
    }
    return 0;
}
#include <errno.h>

// Hulpfunctie om float conversie en validatie te isoleren
static void parse_sd_threshold(const char* arg, double* threshold)
{
    char* endp = NULL;
    errno = 0;
    *threshold = strtod(arg, &endp);
    if (errno == ERANGE || endp == arg || *threshold <= 0.0)
    {
        (void)fprintf(stderr, "invalid float argument for -s '%s'\n", arg);
        {
            exit(-1);
        }
    }
}

// Hulpfunctie voor bestandshandling
static void
handle_file_arg(int opt, FILE** filePtr, const char* arg, const char* mode)
{
    if (checkFileArg(opt, filePtr, arg, mode) != 0)
    {
        exit(-1);
    }
}

static void enforce_uint(int c, unsigned int* target, char* arg)
{
    if (checkUIntArg(c, target, arg) != 0)
    {
        exit(-1);
    }
}

void parse_arguments(int argc, char* argv[], CapConfig* cfg)
{
    int flag;
    while ((flag = getopt(argc, argv, "b:r:z:ht:s:e:c:m:d:w:p:f:D:v:I:lj:")) !=
           -1)
    {
        switch (flag)
        {
        // Group 1: Integrated Numeric Parsers
        case 'j':
            enforce_uint(flag, &cfg->teeth, optarg);
            break;
        case 'f':
            enforce_uint(flag, &cfg->fitN, optarg);
            break;
        case 'e':
            enforce_uint(flag, &cfg->evalue, optarg);
            break;
        case 'v':
            enforce_uint(flag, &cfg->verbose, optarg);
            break;
        case 't':
            enforce_uint(flag, &cfg->time, optarg);
            break;
        case 'z':
            enforce_uint(flag, &cfg->zoom, optarg);
            break;

        // Group 2: Special Logic
        case 'b':
            enforce_uint(flag, &cfg->bph, optarg);
            if (cfg->bph < SECS_HOUR)
            {
                (void)fprintf(stderr, "refusing bph < 3600\n");
                exit(-1);
            }
            break;
        case 'c':
            enforce_uint(flag, &cfg->cvalue, optarg);
            if (cfg->cvalue > CVAL)
            {
                cfg->cvalue = CVAL;
            }
            break;
        case 's':
            parse_sd_threshold(optarg, &cfg->SDthreshold);
            break;
        case 'r':
            if (checkFloatArg(flag, &cfg->rate, optarg) != 0)
            {
                exit(-1);
            }
            break;

        // Group 3: File Handlers
        case 'I':
            handle_file_arg(flag, &cfg->fpInput, optarg, "r");
            break;
        case 'm':
            handle_file_arg(flag, &cfg->fpmaxcor, optarg, "w+");
            break;
        case 'w':
            handle_file_arg(flag, &cfg->fpposition, optarg, "w+");
            break;
        case 'D':
            handle_file_arg(flag, &cfg->fpDefPeak, optarg, "r");
            break;
        case 'p':
            handle_file_arg(flag, &cfg->fptotal, optarg, "w");
            break;

        // Group 4: Flags & Help
        case 'd':
            cfg->device = optarg;
            break;
        case 'l':
            cfg->everyline = EVERY_WIDTH;
            break;
        case 'h':
        default:
            print_usage();
            exit(0);
        }
    }
}
