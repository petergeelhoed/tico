/* 
   A Minimal Capture Program
   This program opens an audio interface for capture, configures it for
   stereo, 16 bit, 44.1kHz, interleaved conventional read/write
   access. Then its reads a chunk of random data from it, and exits. It
   isn't meant to be a real program.
   From on Paul David's tutorial : http://equalarea.com/paul/alsa-audio.html
   Fixes rate and buffer problems
   sudo apt-get install libasound2-dev
   gcc -o alsa-record-example -lasound alsa-record-example.c && ./alsa-record-example hw:0
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <alsa/asoundlib.h>

int main (int argc, char *argv[])
{

    struct timespec beg,end;
    int i;
    int err;
    char *buffer;
    int buffer_frames = 8000;
    unsigned int rate = 48000;
    snd_pcm_t *capture_handle;
    snd_pcm_hw_params_t *hw_params;
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

    if ((err = snd_pcm_open (&capture_handle, argv[1], SND_PCM_STREAM_CAPTURE, 0)) < 0) {
        fprintf (stderr, "cannot open audio device %s (%s)\n", 
                argv[1],
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "audio interface opened\n");

    if ((err = snd_pcm_hw_params_malloc (&hw_params)) < 0) {
        fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params allocated\n");

    if ((err = snd_pcm_hw_params_any (capture_handle, hw_params)) < 0) {
        fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params initialized\n");

    if ((err = snd_pcm_hw_params_set_access (capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
        fprintf (stderr, "cannot set access type (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params access setted\n");

    if ((err = snd_pcm_hw_params_set_format (capture_handle, hw_params, format)) < 0) {
        fprintf (stderr, "cannot set sample format (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params format setted\n");

    if ((err = snd_pcm_hw_params_set_rate_near (capture_handle, hw_params, &rate, 0)) < 0) {
        fprintf (stderr, "cannot set sample rate (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params rate setted\n");

    if ((err = snd_pcm_hw_params_set_channels (capture_handle, hw_params, 1)) < 0) {
        fprintf (stderr, "cannot set channel count (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params channels setted\n");

    if ((err = snd_pcm_hw_params (capture_handle, hw_params)) < 0) {
        fprintf (stderr, "cannot set parameters (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "hw_params setted\n");

    snd_pcm_hw_params_free (hw_params);

    fprintf(stderr, "hw_params freed\n");

    if ((err = snd_pcm_prepare (capture_handle)) < 0) {
        fprintf (stderr, "cannot prepare audio interface for use (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    fprintf(stderr, "audio interface prepared\n");

    buffer = malloc(buffer_frames * snd_pcm_format_width(format) / 8 );

    clock_gettime(CLOCK_REALTIME, &beg);
  //  clock_gettime(CLOCK_REALTIME, &start);
    fprintf(stderr, "buffer allocated %d\n",snd_pcm_format_width(format)); 
    for (i = 0; i < 30*rate/buffer_frames; ++i) {
        if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) {
            fprintf (stderr, "read from audio interface failed %d (%s)\n",
                    err, snd_strerror (err));
            exit (1);
        }
     //   clock_gettime(CLOCK_REALTIME, &end);
        
 //       double f = (double)(end.tv_sec*1e9 + end.tv_nsec);
  //      double g = (double)(start.tv_sec*1e9 + start.tv_nsec);
  //      start = end;
 //       fprintf(stderr,"time %lf\n",(f-g)/1e9 );
         unsigned char lsb;
         signed char msb;
         int in[8000];
         int max = 0;
         int maxpos = 0;
        for (int j = 0; j < buffer_frames*2; j+=2) {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2]=(msb << 8)| lsb ;
            int derivative = (in[j/2]-in[j/2-1] )*(in[j/2]-in[j/2-1] );
            printf("%d\n",derivative );
            if (derivative > max)
            {
            max = derivative;
            maxpos = j/2;
            }
        }

        int width = (maxpos%400)*2/5;
        for (int j = 0; j < width; j++) fprintf(stderr," ");
            fprintf(stderr,"%s\n",i%2==0?"O":"X");
    }

        double f = (double)(end.tv_sec*1e9 + end.tv_nsec);
        double g = (double)(beg.tv_sec*1e9 + beg.tv_nsec);
   //     start = end;
        fprintf(stderr,"Total time %lf\n",(f-g)/1e9 );
    free(buffer);

    fprintf(stderr, "buffer freed\n");

    snd_pcm_close (capture_handle);
    fprintf(stderr, "audio interface closed\n");

    exit (0);
}
