#include "myarr.h"

#include <stdio.h>

void wait(void);
void thread_lock(void);
void thread_unlock(void);
void syncwrite(int* input, unsigned int ArrayLength, char* file);
void* threadWrite(void* inStruct);
void* threadAppend(void* inStruct);
void writearray(int* arr, unsigned int ArrayLength, const char* file);
void printTOD(FILE* out);

unsigned long int syncAppendMyarr(struct myarr* input, FILE* file);
void* threadAppendMyarr(void* inStruct);
