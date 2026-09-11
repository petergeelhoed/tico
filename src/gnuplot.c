/* NOLINTNEXTLINE(cert-dcl37-c,cert-dcl51-cpp) */
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

static void safe_fclose(FILE* filePtr)
{
    if (filePtr && fclose(filePtr))
    {
        perror("fclose");
    }
}

static void print_colored_line(const char* line)
{
    for (const char* p = line; *p; ++p)
    {
        switch (*p)
        {
        case 'X':
            if (fputs(COLOR_RED "X" COLOR_RESET, stdout) == EOF)
            {
                perror("fputs");
            }
            break;

        case 'O':
            if (fputs(COLOR_GREEN "O" COLOR_RESET, stdout) == EOF)
            {
                perror("fputs");
            }
            break;

        default:
            if (fputc(*p, stdout) == EOF)
            {
                perror("fputc");
            }
            break;
        }
    }
}

int gnuplot_cdf(const double* data, size_t length, struct stats* stats)
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

    int retVal = posix_spawnp(&pid, gnuplot_cmd, &actions, NULL, argv, environ);

    posix_spawn_file_actions_destroy(&actions);

    if (retVal != 0)
    {
        // no gnuplot silent continue
        return 0;
    }

    close(stdin_pipe[0]);
    close(stdout_pipe[1]);

    FILE* gnuplot_pipe = fdopen(stdin_pipe[1], "w");

    if (!gnuplot_pipe)
    {
        perror("fdopen");
        //    close(stdout_pipe[0]);
        return 1;
    }
    FILE* gnuplot_out = fdopen(stdout_pipe[0], "r");

    if (!gnuplot_out)
    {
        perror("fdopen");
        safe_fclose(gnuplot_pipe);
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
        const double zval = (data[x] - stats->mean) / stats->stdev;

        printed =
            fprintf(gnuplot_pipe,
                    "%zu %lf %d\n",
                    x,
                    zval,
                    fabs(zval) > STDEV_LIMIT ? pointtype_out : pointtype_in);
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

    safe_fclose(gnuplot_pipe);

    char line[MAX_COLUMNS];

    while (fgets(line, sizeof(line), gnuplot_out))
    {
        print_colored_line(line);
    }

    safe_fclose(gnuplot_out);

    waitpid(pid, NULL, 0);

    return 0;
}
