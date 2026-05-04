
#include <signal.h> // IWYU pragma: export

/**
 * @defgroup mysignal Mysignal Library
 * Signal handling utilities.
 * @{
 */

/** @ingroup mysignal @brief Signal handler for SIGINT (Ctrl+C) and terminal
 * width changes (SIGWINCH)
 *
 * @param signal The signal number received (should be SIGINT)
 */
void sigintHandler(int signal);

/** @ingroup mysignal @brief Set signal handlers for SIGINT and SIGWINCH
 *
 * @details This function sets the signal handlers for SIGINT and SIGWINCH
 * to the sigintHandler function. If there is an error setting the signal
 * handlers, the function will print an error message and exit with an error
 * code.
 */
void setSignalAction(void);

/** @ingroup mysignal @brief Block SIGINT and SIGWINCH signals
 *
 * @param new_set A pointer to a sigset_t structure that will be used to block
 * the SIGINT and SIGWINCH signals. The function will modify this structure to
 * include the blocked signals.
 */
void setupBlockSignals(sigset_t* new_set); // NOLINT(misc-include-cleaner)

/** @} */
