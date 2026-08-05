#include "erf.h"

#include <errno.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#define N 12
#define ISO_LENGTH 40
#define LINESIZE 256
#define HALF 0.5
#define COLOR_RED "\033[31m"
#define COLOR_YELLOW "\033[34m"
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"

extern char** environ;
static int confirm(void)
{
    char line[LINESIZE];

    printf("OK? [Y/n] ");

    if (fgets(line, sizeof(line), stdin) == NULL)
    {
        return 1; /* default yes on EOF */
    }

    return line[0] == '\n' || line[0] == 'y' || line[0] == 'Y';
}

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

static struct stats remove_outliers_and_refit(uint64_t* samples,
                                              unsigned int n,
                                              double limit)
{
    struct stats stats = fit_erf(samples, n);

    const double lower = stats.mean - 2.0 * stats.stdev;
    const double upper = stats.mean + 2.0 * stats.stdev;

    uint64_t filtered[N];
    uint64_t outliers[N];

    unsigned int filtered_n = 0;
    unsigned int outliers_n = 0;

    for (unsigned int i = 0; i < n; i++)
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

    if (outliers_n > 0)
    {
        const double kilo = 1000.0;
        printf("Mean   = %.3f ms\n", stats.mean / kilo);
        printf("Stddev = %s%.3f%s ms\n",
               stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
               stats.stdev / kilo,
               COLOR_RESET);

        for (unsigned int i = 0; i < outliers_n; i++)
        {
            printf("Removed outlier: %lu\n", outliers[i]);
        }
    }

    if (filtered_n > 1)
    {
        return fit_erf(filtered, filtered_n);
    }

    return stats;
}

static double adjust_mean_to_target(double mean, double target)
{
    const double step = (mean < target) ? 5000.0 : -5000.0;

    while (fabs(mean + step - target) < fabs(mean - target))
    {
        mean += step;
    }

    return mean;
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

    printf("\nSorted samples:\n");
    for (int i = 0; i < N; i++)
    {
        printf("%4llu ", (unsigned long long)samples[i]);
    }
    printf("\n");

    const double limit = 150.;
    struct stats stats = remove_outliers_and_refit(samples, N, limit);

    printf("\n\nQQ Gaussian fit:\n");

    stats.mean = adjust_mean_to_target(stats.mean, kilo * value);

    printf("Mean   = %.3f ms\n", stats.mean / kilo);
    printf("Stddev = %s%.3f%s ms\n",
           stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
           stats.stdev / kilo,
           COLOR_RESET);
    printf("prevval = %.3f ms\n", value);

    const double prec_lim = 100.;

    char valuestr[LINESIZE];

    if (stats.stdev < prec_lim)
    {
        int printed =
            snprintf(valuestr, sizeof(valuestr), "%.2f", stats.mean / kilo);
        if (printed < 0 || (size_t)printed >= sizeof(valuestr))
        {
            return EXIT_FAILURE;
        }
    }
    else
    {
        int printed =
            snprintf(valuestr, sizeof(valuestr), "%.1f", stats.mean / kilo);
        if (printed < 0 || (size_t)printed >= sizeof(valuestr))
        {
            return EXIT_FAILURE;
        }
    }

    char arg[LINESIZE];

    if (argc > 1)
    {
        int printed = snprintf(arg, sizeof(arg), "%s # %s", valuestr, argv[1]);
        if (printed < 0 || (size_t)printed >= sizeof(arg))
        {
            return EXIT_FAILURE;
        }
    }
    else
    {
        int printed = snprintf(arg, sizeof(arg), "%s", valuestr);
        if (printed < 0 || (size_t)printed >= sizeof(arg))
        {
            return EXIT_FAILURE;
        }
    }

    char skn_cmd[] = "skn";

    char* argv_exec[] = {skn_cmd, arg, NULL};

    pid_t pid;

    printf("skn '%s'\n", arg);
    if (confirm())
    {
        int retval = posix_spawnp(&pid, "skn", NULL, NULL, argv_exec, environ);

        if (retval == 0)
        {
            int status;
            (void)waitpid(pid, &status, 0);
        }
    }

    return 0;
}
