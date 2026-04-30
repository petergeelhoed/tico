#include <signal.h> // IWYU pragma: export
void sigintHandler(int signal);

void setSignalAction(void);

void setupBlockSignals(sigset_t* new_set); // NOLINT(misc-include-cleaner)

void blockSignal(sigset_t* new_set, sigset_t* old_set);

void unblockSignal(sigset_t* old_set);
