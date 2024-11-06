void syncwrite(int* input, unsigned int NN, char* file);
void syncappend(int* input, unsigned int NN, FILE* file);
void* threadWrite(void* arr);
void* threadAppend(void* arr);
void writearray(int* arr, unsigned int NN, const char* file);
void printTOD(FILE* out);
