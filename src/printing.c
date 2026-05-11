#include "printing.h"

#include "mydefs.h"
#include "mymath.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

extern volatile unsigned int columns;
void printmsg(const char* fmt, ...)
{
    char buf[BUFFER_SIZE];
    va_list args;
    va_start(args, fmt);
    int err = vsnprintf(buf, sizeof(buf), fmt, args);
    if (err < 0)
    {
        strncpy(buf, "[format error]", sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';
    }
    else if ((size_t)err >= sizeof(buf))
    {
        if (sizeof(buf) > 4)
        {
            buf[sizeof(buf) - 4] = '.';
            buf[sizeof(buf) - 3] = '.';
            buf[sizeof(buf) - 2] = '.';
            buf[sizeof(buf) - 1] = '\0';
        }
    }
    va_end(args);

    size_t len = strlen(buf);
    int col = (int)columns - (int)len + 1;
    if (col < 1)
    {
        col = 1;
    }
    (void)fprintf(stderr,
                  "\033[s\033[12;0H\033[%dG\033[96m%s\033[0m\033[u",
                  col,
                  buf);
}

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds)
{
    if (everyline)
    {
        static char tmp[DOUBLE_BUF];
        (void)sprintf(tmp, "%4.2f", beatError);
        (void)sprintf(tmp + BEAT_WIDTH - 1, "ms%+5.1f", fittedRate);
        (void)sprintf(tmp + BEAT_WIDTH + RATE_WIDTH - 2, "s/d");
        tmp[EVERY_WIDTH] = '\0';
        (void)fprintf(stderr, "%s", tmp);
    }
    else
    {
        (void)fprintf(
            stderr,
            "\033[s\033[2;0H\033[0K%8.2fms   %9.1fs/d   %12.2fs\033[u",
            beatError,
            fittedRate,
            seconds);
    }
}

void printspaces(int maxpos,
                 size_t hexvalue,
                 size_t mod,
                 size_t acolumns,
                 double avgPos,
                 size_t correlationThreshold)
{
    while (maxpos < (int)mod)
    {
        maxpos += (int)mod;
    }
    while (avgPos < (double)mod)
    {
        avgPos += (double)mod;
    }
    acolumns = acolumns > MAX_COLUMNS ? DEFAULT_COLUMNS : acolumns;
    size_t width = modSigned(maxpos, mod) * acolumns / mod;
    size_t widtha = modSigned((int)lround(avgPos), mod) * acolumns / mod;

    char spaces[MAX_COLUMNS];
    memset(spaces, ' ', width);
    spaces[width] = '\0';
    if (widtha < width)
    {
        spaces[widtha] = '|';
    }

    (void)fprintf(stderr,
                  "%s%s%zX\033[0m",
                  spaces,
                  hexvalue < correlationThreshold ? "\033[31m" : "\033[32m",
                  hexvalue);

    memset(spaces, ' ', width);
    if (widtha > width)
    {
        spaces[widtha - width - 1] = '|';
        spaces[widtha - width] = '\0';
        (void)fprintf(stderr, "%s", spaces);
    }
    (void)fprintf(stderr, "\n");
}
