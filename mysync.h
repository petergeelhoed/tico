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
