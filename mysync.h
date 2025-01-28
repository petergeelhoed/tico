#include "myarr.h"

void wait();
void thread_lock();
void thread_unlock();
void syncwrite(int* input, unsigned int NN, char* file);
long unsigned int syncappend(int* input, unsigned int NN, FILE* file);
void* threadWrite(void* arr);
void* threadAppend(void* arr);
void writearray(int* arr, unsigned int NN, const char* file);
void printTOD(FILE* out);

void syncappendDouble(double* input, unsigned int NN, FILE* file);
void* threadAppendDouble(void* arr);

unsigned long int syncAppendMyarr(struct myarr* input, FILE* file);
void* threadAppendMyarr(void* arr);
