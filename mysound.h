#include <alsa/asoundlib.h>
#include <fftw3.h>

snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate);
int readBuffer(snd_pcm_t* capture_handle,
               int NN,
               char* buffer,
               int* derivative);
void readBufferRaw(snd_pcm_t* capture_handle, int NN, char* buffer, int* in);

int readShiftedBuffer(int* derivative,
                      snd_pcm_t* capture_handle,
                      int NN,
                      char* buffer,
                      int maxpos,
                      int* totalshift, FILE* fpInput);
