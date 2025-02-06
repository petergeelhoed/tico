#include "parseargs.h"
#include "mylib.h"
#include <unistd.h>

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
            if (checkUIntArg(c, bph, optarg) != 0 || *bph < 3600)
            {
                fprintf(stderr, "refusing bph < 3600 %d\n", *bph);
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
