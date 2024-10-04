void syncwrite(int* input, int NN, char* file);
void syncappend(int* input, int NN, FILE* file);
void* threadWrite(void* arr);
void* threadAppend(void* arr);
void writearray(int* arr, int NN, const char* file);
