#include "myarr.h"
#include <alsa/asoundlib.h>
#include <fftw3.h>

snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate);
int readBuffer(snd_pcm_t* capture_handle,
               unsigned int NN,
               char* buffer,
               int* derivative);
void readBufferRaw(snd_pcm_t* capture_handle,
                   unsigned int NN,
                   char* buffer,
                   int* in);

int readBufferOrFile(int* derivative,
                     snd_pcm_t* capture_handle,
                     unsigned int NN,
                     char* buffer,
                     FILE* fpInput);
int getData(FILE* rawfile,
            FILE* fpInput,
            snd_pcm_t* capture_handle,
            snd_pcm_format_t format,
            char* device,
            unsigned int rate,
            char* buffer,
            struct myarr derivative);
