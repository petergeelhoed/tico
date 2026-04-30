#include <signal.h> // IWYU pragma: export
void sigintHandler(int signal);

void setSignalAction(void);

void setupBlockSignals(sigset_t* new_set); // NOLINT(misc-include-cleaner)
