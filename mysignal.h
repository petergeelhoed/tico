#include <signal.h>
void sigint_handler(int signal);

void set_signal_action(void);

void setup_block_signals(sigset_t* new_set);

void block_signal(sigset_t* new_set, sigset_t* old_set);

void unblock_signal(sigset_t* old_set);
