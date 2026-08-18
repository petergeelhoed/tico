#include "erf.h"

#include <errno.h>
#include <math.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#define N 12
#define ISO_LENGTH 40
#define LINESIZE 512
#define VALUESIZE 512
#define COLOR_RED "\033[31m"
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"

extern char** environ;

static int compare_double(const void* first, const void* second)
{
    double dfirst = *(const double*)first;
    double dsecond = *(const double*)second;

    return (dfirst > dsecond) - (dfirst < dsecond);
}

static int is_effectively_zero(double value)
{
    const double epsilon = 1e-9;
    return fabs(value) < epsilon;
}

static int confirm(double* adjust)
{
    char line[LINESIZE];

    *adjust = 0.0;

    printf("OK? [Y/n or seconds to add] ");

    if (fgets(line, sizeof(line), stdin) == NULL)
    {
        return 1;
    }

    if (line[0] == '\n' || line[0] == 'y' || line[0] == 'Y')
    {
        return 1;
    }

    char* end;
    errno = 0;
    double value = strtod(line, &end);

    if ((end != line) && (*end == '\n' || *end == '\0') && (errno != ERANGE))
    {
        *adjust = value;
    }

    return 0;
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
        exit(EXIT_FAILURE);
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

static void unwrap(double* data, size_t n)
{
    const double modulo = 5.0;
    const double half = modulo / 2.0;

    for (size_t i = 0; i < n - 1; ++i)
    {
        if (data[i + 1] - data[i] > half)
        {
            /* wrap is between i and i+1 */
            for (size_t j = 0; j <= i; ++j)
            {
                data[j] += modulo;
            }
            break;
        }
    }
}

static struct stats remove_outliers_and_refit(double* samples,
                                              unsigned int n,
                                              double limit)
{
    struct stats stats = fit_erf(samples, n);

    const double lower = stats.mean - 1.5 * stats.stdev;
    const double upper = stats.mean + 1.5 * stats.stdev;

    double filtered[N];
    double outliers[N];

    unsigned int filtered_n = 0;
    unsigned int outliers_n = 0;

    for (unsigned int i = 0; i < n; i++)
    {
        if (samples[i] >= lower && samples[i] <= upper)
        {
            filtered[filtered_n++] = samples[i];
            if (!is_effectively_zero(stats.stdev))
            {
                printf("%7.3f %7.3f\n",
                       samples[i],
                       (stats.mean - samples[i]) / stats.stdev);
            }
            else
            {
                printf("%7.3f %7.3f\n", samples[i], 0.0);
            }
        }
        else
        {
            outliers[outliers_n++] = samples[i];
            if (!is_effectively_zero(stats.stdev))
            {
                printf("%7.3f" COLOR_RED " %7.3f" COLOR_RESET "\n",
                       samples[i],
                       (stats.mean - samples[i]) / stats.stdev);
            }
            else
            {
                printf("%7.3f" COLOR_RED " %7.3f" COLOR_RESET "\n",
                       samples[i],
                       0.0);
            }
        }
    }

    if (outliers_n > 0)
    {
        printf("Mean   = %.3f s\n", stats.mean);
        printf("Stddev = %s%.3f%s s\n",
               stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
               stats.stdev,
               COLOR_RESET);

        for (unsigned int i = 0; i < outliers_n; i++)
        {
            if (!is_effectively_zero(stats.stdev))
            {
                printf("Removed outlier: %7.3f %7.3f\n",
                       outliers[i],
                       (stats.mean - outliers[i]) / stats.stdev);
            }
            else
            {
                printf("Removed outlier: %7.3f %7.3f\n", outliers[i], 0.0);
            }
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
    const double step = (mean < target) ? 5.0 : -5.0;

    while (fabs(mean + step - target) < fabs(mean - target))
    {
        mean += step;
        printf("%.3f %.3f\n", mean, target);
    }

    return mean;
}

static void build_arg(char* valuestr,
                      char* arg,
                      int argc,
                      char** argv,
                      struct stats stats)
{
    const double prec_lim = .100;
    if (stats.stdev < prec_lim)
    {
        int printed = snprintf(valuestr, VALUESIZE, "%.2f", stats.mean);
        if (printed < 0 || (size_t)printed >= VALUESIZE)
        {
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        int printed = snprintf(valuestr, VALUESIZE, "%.1f", stats.mean);
        if (printed < 0 || (size_t)printed >= VALUESIZE)
        {
            exit(EXIT_FAILURE);
        }
    }

    if (argc > 1)
    {
        int printed = snprintf(arg, LINESIZE, "%s # %s", valuestr, argv[1]);
        if (printed < 0 || (size_t)printed >= LINESIZE)
        {
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        int printed = snprintf(arg, LINESIZE, "%s", valuestr);
        if (printed < 0 || (size_t)printed >= LINESIZE)
        {
            exit(EXIT_FAILURE);
        }
    }
}

static void get_values(double* samples)
{
    struct timespec tspec;

    printf("Press ENTER %d times, every 5 seconds\n\n", N);

    for (int i = 0; i < N; i++)
    {
        printf("[%d/%d] ", i + 1, N);

        while (getchar() != '\n')
        {
        }

        int err = clock_gettime(CLOCK_REALTIME, &tspec);
        if (err)
        {
            exit(EXIT_FAILURE);
        }
        const long long mod = 5;
        const double NANO = 1e-9;
        samples[i] =
            -(double)(tspec.tv_sec % mod) - (double)tspec.tv_nsec * NANO;

        struct tm tm_info;
        char iso8601[ISO_LENGTH];
        if (NULL == gmtime_r(&tspec.tv_sec, &tm_info))
        {
            exit(EXIT_FAILURE);
        }
        if (0 ==
            strftime(iso8601, sizeof(iso8601), "%Y-%m-%dT%H:%M:%S", &tm_info))
        {
            exit(EXIT_FAILURE);
        }

        const long MEGA = 1000000L;
        printf("%.3f s  %s.%03ldZ\n",
               samples[i],
               iso8601,
               tspec.tv_nsec / MEGA);
    }
}

int main(int argc, char** argv)
{
    double value = read_last_value("/var/www/temp/seiko");

    double samples[N];
    get_values(samples);

    const double limit = .150;
    qsort(samples, N, sizeof(double), compare_double);
    unwrap(samples, N);

    struct stats stats = remove_outliers_and_refit(samples, N, limit);
    printf("\n\nQQ Gaussian fit:\n");

    stats.mean = adjust_mean_to_target(stats.mean, value);

    printf("Mean   = %.3f s\n", stats.mean);
    printf("Stddev = %s%.3f%s s\n",
           stats.stdev > limit ? COLOR_RED : COLOR_GREEN,
           stats.stdev,
           COLOR_RESET);
    printf("prevval = %.3f s\n", value);

    char valuestr[VALUESIZE];
    char arg[LINESIZE];

    build_arg(valuestr, arg, argc, argv, stats);
    char skn_cmd[] = "skn";
    char* argv_exec[] = {skn_cmd, arg, NULL};

    pid_t pid = -1;

    for (;;)
    {
        printf("skn '%s'\n", arg);

        double adjust = 0.0;

        if (confirm(&adjust))
        {
            break;
        }

        if (fabs(adjust) > 0.0)
        {
            stats.mean += adjust;
            build_arg(valuestr, arg, argc, argv, stats);
            continue;
        }

        return EXIT_FAILURE;
    }

    int retval = posix_spawnp(&pid, "skn", NULL, NULL, argv_exec, environ);

    if (retval == 0)
    {
        int status;
        retval = waitpid(pid, &status, 0);
    }

    return retval;
}
