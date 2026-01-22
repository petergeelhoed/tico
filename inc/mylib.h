#pragma once
#include "config.h"
#include "crosscorint.h"
#include "myarr.h"
#include "resources.h"
#include <stdio.h>

void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds);

void printspaces(int maxpos,
                 double hexvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double avg_pos,
                 unsigned int correlationThreshold);

void writefile(FILE* filePtr, int* array, unsigned int ArrayLength);
void writefileDouble(FILE* filePtr, double* array, unsigned int ArrayLength);

double getBeatError(const struct myarr* totalTick, double rate, int verbose);

void fillReference(FILE* fpDefPeak,
                   struct myarr* reference,
                   unsigned int teeth);

void print_finals(CapConfig* cfg,
                  AppResources* res,
                  unsigned int ArrayLength,
                  unsigned int totalTickTock,
                  int toothshift);
void shift_buffer_data(unsigned int* ticktock,
                       struct myarr* subpos,
                       struct myarr* maxpos,
                       struct myarr* maxvals);
