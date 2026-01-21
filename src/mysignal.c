#include "mysignal.h"

#include "mydefs.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

extern int keepRunning;
extern unsigned int columns;

void sigint_handler(int signal)
{
    if (signal == SIGINT)
    {
        keepRunning = 0;
    }
    else if (signal == SIGWINCH)
    {
        struct winsize windowSize;
        ioctl(STDOUT_FILENO,
              TIOCGWINSZ,
              &windowSize); // NOLINT(misc-include-cleaner)
        columns = (unsigned int)windowSize.ws_col;
        (void)fprintf(stderr, "new width %d\n", columns);
    }
    else
    {
        (void)raise(signal);
    }
}

void set_signal_action(void)
{
    struct sigaction sact;
    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = sigint_handler;

    // Set signal handlers
    if (sigaction(SIGWINCH, &sact, NULL) == -1)
    {
        perror("failed sigaction SIGWINCH");
        exit(ERROR_SIGNAL);
    }
    if (sigaction(SIGINT, &sact, NULL) == -1)
    {
        perror("failed sigaction SIGINT");
        exit(ERROR_SIGNAL);
    }
}

void setup_block_signals(sigset_t* new_set) // NOLINT(misc-include-cleaner)
{
    struct sigaction sact;

    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = sigint_handler;
    if (sigaction(SIGWINCH, &sact, NULL) != 0)
    {
        (void)fprintf(stderr, "sigaction() error");
        exit(ERROR_SIGNAL);
    }
    sigemptyset(new_set);
    sigaddset(new_set, SIGWINCH);
    sigaddset(new_set, SIGINT);
}

void block_signal(sigset_t* new_set, sigset_t* old_set)
{
    if (sigprocmask(SIG_BLOCK, new_set, old_set) != 0)
    {
        (void)fprintf(stderr, "block sigprocmask() error");
        exit(ERROR_SIGNAL);
    }
}

void unblock_signal(sigset_t* old_set)
{
    if (sigprocmask(SIG_SETMASK, old_set, NULL) != 0)
    {
        (void)fprintf(stderr, "unblock sigprocmask() error");
        exit(ERROR_SIGNAL);
    }
}
