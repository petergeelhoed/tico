#pragma once
#include "stats.h"

#include <unistd.h>

int gnuplot_cdf(const double* data, size_t length, struct stats* stats);
