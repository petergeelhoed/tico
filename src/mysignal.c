
#include "mysignal.h"

#include "appstate.h"
#include "mydefs.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

static struct AppState* s_appState = NULL;

void sigintHandler(int signal)
{
    if (!s_appState)
    {
        return;
    }
    if (signal == SIGINT)
    {
        s_appState->keepRunning = 0;
    }
    else if (signal == SIGWINCH)
    {
        struct winsize windowSize;
        ioctl(STDOUT_FILENO,
              TIOCGWINSZ,
              &windowSize); // NOLINT(misc-include-cleaner)
        s_appState->columns = (unsigned int)windowSize.ws_col;
        (void)fprintf(stderr, "new width %d\n", s_appState->columns);
    }
    else
    {
        (void)raise(signal);
    }
}

void setSignalAction(struct AppState* appState)
{
    s_appState = appState;
    struct sigaction sact;
    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = sigintHandler;

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

void setupBlockSignals(sigset_t* new_set) // NOLINT(misc-include-cleaner)
{
    struct sigaction sact;

    sigemptyset(&sact.sa_mask);
    sact.sa_flags = 0;
    sact.sa_handler = sigintHandler;
    if (sigaction(SIGWINCH, &sact, NULL) != 0)
    {
        (void)fprintf(stderr, "sigaction() error");
        exit(ERROR_SIGNAL);
    }
    sigemptyset(new_set);
    sigaddset(new_set, SIGWINCH);
    sigaddset(new_set, SIGINT);
}
