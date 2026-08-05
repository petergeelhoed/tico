#include "erf.h"

#include <stdio.h>
#include <time.h>

#define N 12
#define ISO_LENGTH 40

int main(void)
{
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
    printf("\n\nQQ Gaussian fit:\n");
    printf("Mean   = %.3f ms\n", stats.mean);
    printf("Stddev = %.3f ms\n", stats.stdev);

    return 0;
}
