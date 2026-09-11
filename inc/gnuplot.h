#pragma once
#include "stats.h"

#include <unistd.h>

int gnuplot_cdf(double* data, size_t length, struct stats* stats);
