#include <alsa/asoundlib.h>
#include <fftw3.h>
#include "myarr.h"

snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate);
int readBuffer(snd_pcm_t* capture_handle,
               unsigned int NN,
               char* buffer,
               int* derivative);
void readBufferRaw(snd_pcm_t* capture_handle,
                   unsigned int NN,
                   char* buffer,
                   int* in);

int readShiftedBuffer(int* derivative,
                      snd_pcm_t* capture_handle,
                      unsigned int NN,
                      char* buffer,
                      int maxpos,
                      FILE* fpInput);
int getData(unsigned int maxp,
            int* totalshift,
            FILE* rawfile,
            FILE* fpInput,
            snd_pcm_t* capture_handle,
            snd_pcm_format_t format,
            char* device,
            unsigned int rate,
            char* buffer,
            struct myarr derivative,
            unsigned int totalI);
