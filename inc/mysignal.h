#include <signal.h> // IWYU pragma: export
void sigint_handler(int signal);

void set_signal_action(void);

void setup_block_signals(sigset_t* new_set); // NOLINT(misc-include-cleaner)

void block_signal(sigset_t* new_set, sigset_t* old_set);

void unblock_signal(sigset_t* old_set);
