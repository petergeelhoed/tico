#include "printing.h"

#include "modsigned.h"
#include "mydefs.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

// Terminal color and control escape sequences
#define COLOR_CYAN "\033[96m"
#define COLOR_RED "\033[31m"
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"
#define CURSOR_SAVE "\033[s"
#define CURSOR_RESTORE "\033[u"
#define CLEAR_LINE "\033[0K"
#define CURSOR_ROW12 "\033[12;0H"

#define CURSOR_ROW2 "\033[2;0H"
#define CURSOR_COL_FMT "\033[%dG"

extern volatile unsigned int columns;
void printmsg(const char* fmt, ...)
{
    char buf[BUFFER_SIZE];
    va_list args;
    va_start(args, fmt);
    // NOLINTBEGIN(clang-analyzer-valist.Uninitialized)
    int err = vsnprintf(buf, sizeof(buf), fmt, args);
    // NOLINTEND(clang-analyzer-valist.Uninitialized)
    va_end(args);

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

    size_t len = strlen(buf);
    int col = (int)columns - (int)len + 1;
    if (col < 1)
    {
        col = 1;
    }
    print(CURSOR_SAVE CURSOR_ROW12 CURSOR_COL_FMT COLOR_CYAN
          "%s" COLOR_RESET CURSOR_RESTORE,
          col,
          buf);
}

void print(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    // NOLINTBEGIN(clang-analyzer-valist.Uninitialized)
    (void)vfprintf(stderr, fmt, args);
    // NOLINTEND(clang-analyzer-valist.Uninitialized)
    va_end(args);
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
        print("%s", tmp);
    }
    else
    {
        print(CURSOR_SAVE CURSOR_ROW2 CLEAR_LINE
              "%8.2fms   %9.1fs/d   %12.2fs" CURSOR_RESTORE,
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

    print("%s%s%zX" COLOR_RESET,
          spaces,
          hexvalue < correlationThreshold ? COLOR_RED : COLOR_GREEN,
          hexvalue);

    memset(spaces, ' ', width);
    if (widtha > width)
    {
        spaces[widtha - width - 1] = '|';
        spaces[widtha - width] = '\0';
        print("%s", spaces);
    }
    print("\n");
}
