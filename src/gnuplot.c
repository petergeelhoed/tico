#define _POSIX_C_SOURCE 200809L

#include "gnuplot.h"

#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

int gnuplot_cdf(double* data, size_t length)
{
    int filedescriptors[2];

    if (pipe(filedescriptors) != 0)
    {
        perror("pipe");
        return 1;
    }

    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);

    /* Child stdin <- pipe read end */
    posix_spawn_file_actions_adddup2(&actions,
                                     filedescriptors[0],
                                     STDIN_FILENO);

    /* Close unneeded fds in child */
    posix_spawn_file_actions_addclose(&actions, filedescriptors[0]);
    posix_spawn_file_actions_addclose(&actions, filedescriptors[1]);

    pid_t pid;
    char gnuplot_cmd[] = "gnuplot";
    char* argv[] = {gnuplot_cmd, NULL};

    if (posix_spawnp(&pid, gnuplot_cmd, &actions, NULL, argv, environ) != 0)
    {
        perror("posix_spawnp");
        return 1;
    }

    posix_spawn_file_actions_destroy(&actions);

    /* Parent only writes */
    close(filedescriptors[0]);

    FILE* gnuplot_pipe = fdopen(filedescriptors[1], "w");
    if (!gnuplot_pipe)
    {
        perror("fdopen");
        return 1;
    }

    int printed =
        fprintf(gnuplot_pipe,
                //"set term dumb; uns key; uns xtics; uns ytics; plot "
                "set term dumb; uns key; plot "
                "'-' with points pt 15\n");
    if (printed < 0)
    {
        perror("pipe");
    }

    for (size_t x = 0; x < length; ++x)
    {
        printed = fprintf(gnuplot_pipe, "%lf\n", *(++data));
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

    waitpid(pid, NULL, 0);

    return 0;
}
