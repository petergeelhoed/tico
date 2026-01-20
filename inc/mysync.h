#include "myarr.h"

#include <stdio.h>

void wait(void);
int syncwrite(int* input, unsigned int ArrayLength, const char* file);
void writearray(int* arr, unsigned int ArrayLength, const char* file);
void printTOD(FILE* out);

void syncAppendMyarr(struct myarr* input, FILE* file);
void* threadAppendMyarr(void* inStruct);
