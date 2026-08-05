#include "erf.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 12
#define ISO_LENGTH 40
#define LINESIZE 256
#define HALF 0.5
#define COLOR_RED "\033[31m"
#define COLOR_YELLOW "\033[34m"
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"

static double read_last_value(const char* filename)
{
    FILE* fileptr = fopen(filename, "r");
    if (fileptr == NULL)
    {
        perror(filename);
        exit(EXIT_FAILURE);
    }

    char line[LINESIZE];
    char last[LINESIZE] = "";

    while (fgets(line, sizeof(line), fileptr) != NULL)
    {
        if (line[0] != '\n')
        {
            strncpy(last, line, sizeof(last) - 1);
            last[sizeof(last) - 1] = '\0';
        }
    }

    if (fclose(fileptr))
    {
        exit(-1);
    }

    char timestamp[LINESIZE];
    char value_str[LINESIZE];

    if (sscanf(last, "%63s %63s", timestamp, value_str) != 2)
    {
        exit(EXIT_FAILURE);
    }

    char* endptr = NULL;
    errno = 0;
    double value = strtod(value_str, &endptr);

    if ((endptr == value_str) || (*endptr != '\0') || (errno == ERANGE))
    {
        exit(EXIT_FAILURE);
    }

    return value;
}

int main(int argc, char** argv)
{
    double value = read_last_value("/var/www/temp/seiko");

    uint64_t samples[N];
    struct timespec tspec;
    const uint64_t mod = 5000;
    const uint64_t kilo = 1000ULL;
    const uint64_t mega = 1000000ULL;
    const long MEGA = 1000000UL;
    printf("Press ENTER %d times\n\n", N);

    for (int i = 0; i < N; i++)
    {
        printf("[%d/%d] ", i + 1, N);

        while (getchar() != '\n')
        {
        }

        int err = clock_gettime(CLOCK_REALTIME, &tspec);
        if (err)
        {
            exit(-1);
        }

        uint64_t msec =
            (uint64_t)tspec.tv_sec * kilo + (uint64_t)tspec.tv_nsec / mega;

        samples[i] = msec % mod;
        struct tm time;
        char iso8601[ISO_LENGTH];
        if (NULL == gmtime_r(&tspec.tv_sec, &time))
        {
            exit(-1);
        }
        if (0 == strftime(iso8601, sizeof(iso8601), "%Y-%m-%dT%H:%M:%S", &time))
        {
            exit(-1);
        }
        printf("%llu ms  %s.%03ldZ\n",
               (unsigned long long)samples[i],
               iso8601,
               tspec.tv_nsec / MEGA);
    }

    struct stats stats = fit_erf(samples, N);
    printf("\nSorted samples:\n");
    for (int i = 0; i < N; i++)
    {
        printf("%4llu ", (unsigned long long)samples[i]);
    }
    printf("\n");

    const double limit = 150.;
    printf("Mean   = %.3f ms\n", stats.mean / kilo);
    printf("Stddev = %s%.3f%s ms\n",
           stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
           stats.stdev / kilo,
           COLOR_RESET);

    /* Remove points outside ±2σ */
    const double lower = stats.mean - 2.0 * stats.stdev;
    const double upper = stats.mean + 2.0 * stats.stdev;

    uint64_t filtered[N];
    uint64_t outliers[N];
    unsigned int filtered_n = 0;
    unsigned int outliers_n = 0;

    for (int i = 0; i < N; i++)
    {
        if ((double)samples[i] >= lower && (double)samples[i] <= upper)
        {
            filtered[filtered_n++] = samples[i];
        }
        else
        {
            outliers[outliers_n++] = samples[i];
        }
    }

    for (unsigned int i = 0; i < outliers_n; i++)
    {
        printf("Removed outlier: %lu\n", outliers[i]);
    }
    /* Refit using filtered data */
    if (filtered_n > 1)
    {
        stats = fit_erf(filtered, filtered_n);
    }
    printf("\n\nQQ Gaussian fit:\n");

    const double target = kilo * value;
    const double step = (stats.mean < target) ? 5000.0 : -5000.0;
    while (fabs(stats.mean + step - target) < fabs(stats.mean - target))
    {
        stats.mean += step;
    }

    printf("Mean   = %.3f ms\n", stats.mean / kilo);
    printf("Stddev = %s%.3f%s ms\n",
           stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
           stats.stdev / kilo,
           COLOR_RESET);
    printf("prevval = %.3f ms\n", value);

    const double prec_lim = 100.;
    char extra[LINESIZE] = "";
    if (argc > 1)
    {
        extra[0] = '#';
        strncpy(extra + 1, argv[1], LINESIZE - 2);
        extra[LINESIZE - 1] = '\0';
    }

    if (stats.stdev < prec_lim)
    {
        printf("\nskn %.2f %s\n", stats.mean / kilo, extra);
    }
    else
    {
        printf("\nskn %.1f %s\n", stats.mean / kilo, extra);
    }
    return 0;
}
