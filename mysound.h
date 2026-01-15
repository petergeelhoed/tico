#include "myarr.h"
#include <alsa/asoundlib.h> // IWYU pragma: export
#include <stdio.h>

// NOLINTNEXTLINE(misc-include-cleaner)
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate);
int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative);
void readBufferRaw(snd_pcm_t* capture_handle,
                   char* buffer,
                   struct myarr* data_in);

int readBufferOrFile(int* derivative,
                     snd_pcm_t* capture_handle,
                     unsigned int ArrayLength,
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
