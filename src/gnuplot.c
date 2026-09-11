#define _POSIX_C_SOURCE 200809L

#include "gnuplot.h"

#include "mydefs.h"
#include "stats.h"

#include <math.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

int gnuplot_cdf(double* data, size_t length, struct stats* stats)
{

    struct winsize windowSize;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &windowSize);
    int columns = windowSize.ws_col;

    int stdin_pipe[2];
    int stdout_pipe[2];

    if (pipe(stdin_pipe) != 0)
    {
        perror("pipe");
        return 1;
    }

    if (pipe(stdout_pipe) != 0)
    {
        perror("pipe");
        close(stdin_pipe[0]);
        close(stdin_pipe[1]);
        return 1;
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);

    /* child stdin <- parent */
    posix_spawn_file_actions_adddup2(&actions, stdin_pipe[0], STDIN_FILENO);

    /* child stdout -> parent */
    posix_spawn_file_actions_adddup2(&actions, stdout_pipe[1], STDOUT_FILENO);

    /* capture stderr too */
    posix_spawn_file_actions_adddup2(&actions, stdout_pipe[1], STDERR_FILENO);

    posix_spawn_file_actions_addclose(&actions, stdin_pipe[1]);
    posix_spawn_file_actions_addclose(&actions, stdout_pipe[0]);

    pid_t pid;
    char gnuplot_cmd[] = "gnuplot";
    char* argv[] = {gnuplot_cmd, NULL};

    int rc = posix_spawnp(&pid, gnuplot_cmd, &actions, NULL, argv, environ);

    posix_spawn_file_actions_destroy(&actions);

    if (rc != 0)
    {
        fprintf(stderr, "Failed to start gnuplot: %s\n", strerror(rc));
        return 1;
    }

    close(stdin_pipe[0]);
    close(stdout_pipe[1]);

    FILE* gnuplot_pipe = fdopen(stdin_pipe[1], "w");
    FILE* gnuplot_out = fdopen(stdout_pipe[0], "r");

    if (!gnuplot_pipe)
    {
        perror("fdopen");
        close(stdout_pipe[0]);
        return 1;
    }

    if (!gnuplot_out)
    {
        perror("fdopen");
        fclose(gnuplot_pipe);
        return 1;
    }

    int printed = fprintf(
        gnuplot_pipe,
        //"set term dumb; uns key; uns xtics; uns ytics; plot "
        "set term dumb size %d,24; uns key; unset xtics; set ytics 1 out; plot "
        "[-1:%lu][%lf:%lf]"
        "'-' u 1:2:3 with points pt var\n",
        columns,
        length,
        (data[0] - (stats->mean)) / stats->stdev - 1,
        (data[length - 1] - (stats->mean)) / stats->stdev + 1);
    if (printed < 0)
    {
        perror("pipe");
    }

    const int pointtype_out = 24;
    const int pointtype_in = 15;

    for (size_t x = 0; x < length; ++x)
    {
        double z = (data[x] - stats->mean) / stats->stdev;

        printed = fprintf(gnuplot_pipe,
                          "%zu %lf %d\n",
                          x,
                          z,
                          fabs(z) > STDEV_LIMIT ? pointtype_out : pointtype_in);
        if (printed < 0)
        {
            perror("pipe");
        }
    }

    printed = fprintf(gnuplot_pipe, "e\n");
    if (printed < 0)
    {
        perror("pipe");
    }
    if (fflush(gnuplot_pipe))
    {
        perror("flush");
    }

    if (fclose(gnuplot_pipe))
    {
        perror("fclose");
    }

    char line[MAX_COLUMNS];

    while (fgets(line, sizeof(line), gnuplot_out))
    {
        for (char* p = line; *p; ++p)
        {
            if (*p == 'X')
            {
                if (EOF == fputs(COLOR_RED "X" COLOR_RESET, stdout))
                {
                    perror("fputs");
                }
            }
            else if (*p == 'O')
            {
                if (EOF == fputs(COLOR_GREEN "O" COLOR_RESET, stdout))
                {
                    perror("fputs");
                }
            }
            else
            {
                if (EOF == fputc(*p, stdout))
                {
                    perror("fputc");
                }
            }
        }
    }

    if (fclose(gnuplot_out))
    {
        perror("fclose");
    }

    waitpid(pid, NULL, 0);

    return 0;
}
