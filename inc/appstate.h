#pragma once

/** @brief Holds application-wide runtime state. */
typedef struct AppState
{
    /** Nonzero while the application should continue running. */
    volatile int keepRunning;
    /** Current terminal width in columns. */
    volatile unsigned int columns;
} AppState;
