#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>

int main (int argc, char *argv[])
{
    unsigned int rate = 48000;
    int bph = 21600;
    int buffer_frames = rate*3600/bph;
    int mvalue = 20;
    int time = 30;
    int c;
    int verbose = 0;
    double threshold =3.;
    while ((c = getopt (argc, argv, "b:r:m:ht:vs:")) != -1)
    {
        switch (c)
        {
            case 'v':
                verbose = 1;
                break;
            case 's':
                threshold = atoi(optarg);
                break;
            case 't':
                time = atoi(optarg);
                break;
            case 'b':
                bph = atoi(optarg);
                break;
            case 'm':
                mvalue = atoi(optarg);
                break;
            case 'r':
                rate = atoi(optarg);
                break;
            case 'h':
                fprintf (stderr, "usage:\n capture device (default default:1)\noptions:\n -m <fraction of tick to modulate and plot (default: 20)\n -b bph of the watch (default: 21600/h) \n -r sampling rate (default: 48000Hz)\n -t <measurment time> (default: 30s)\n -v verbose, print points to stdout\n time, tick position modulo(3600/rate), deviation, σ\n -s cutoff standarddeviation (default: 3)\n"); 
                exit(0);

            default:
                break;
        }
    }
    int i;
    int err;
    char *buffer;
    int maxes[time * rate/buffer_frames];
    int mod = buffer_frames/mvalue;

    snd_pcm_t *capture_handle;
    snd_pcm_hw_params_t *hw_params;
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;

    if ((err = snd_pcm_open (&capture_handle, argv[optind], SND_PCM_STREAM_CAPTURE, 0)) < 0) {
        fprintf (stderr, "cannot open audio device %s (%s)\n", 
                argv[1],
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_malloc (&hw_params)) < 0) {
        fprintf (stderr, "cannot allocate hardware parameter structure (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_any (capture_handle, hw_params)) < 0) {
        fprintf (stderr, "cannot initialize hardware parameter structure (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_set_access (capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
        fprintf (stderr, "cannot set access type (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_set_format (capture_handle, hw_params, format)) < 0) {
        fprintf (stderr, "cannot set sample format (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_set_rate_near (capture_handle, hw_params, &rate, 0)) < 0) {
        fprintf (stderr, "cannot set sample rate (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params_set_channels (capture_handle, hw_params, 1)) < 0) {
        fprintf (stderr, "cannot set channel count (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    if ((err = snd_pcm_hw_params (capture_handle, hw_params)) < 0) {
        fprintf (stderr, "cannot set parameters (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    snd_pcm_hw_params_free (hw_params);

    if ((err = snd_pcm_prepare (capture_handle)) < 0) {
        fprintf (stderr, "cannot prepare audio interface for use (%s)\n",
                snd_strerror (err));
        exit (1);
    }

    buffer = malloc(buffer_frames * snd_pcm_format_width(format) / 8 );

    FILE *fp;
    char out[16];
    int wdth;
    fp =popen("/usr/bin/tput cols" , "r");
    fgets(out,16,fp);
    fclose(fp);
    wdth=atoi(out);

    fprintf(stderr, "Found COLUMNS=%d, width = %.3fms  /  %.1fμs/character\n",
            wdth - 1,
            mod*1000./rate,
            mod*1000000./rate/(wdth-1));

    for (i = 0; i < time * rate/buffer_frames; ++i)
    {
        if ((err = snd_pcm_readi (capture_handle, buffer, buffer_frames)) != buffer_frames) {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
        unsigned char lsb;
        signed char msb;
        int in[buffer_frames];
        int max = 0;
        int maxpos = 0;
        for (int j = 0; j < buffer_frames*2; j+=2) {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2]=(msb << 8)| lsb ;
            int derivative = (in[j/2]-in[j/2-1] )*(in[j/2]-in[j/2-1] );
            //  printf("%d\n",derivative );
            if (derivative > max)
            {
                max = derivative;
                maxpos = j/2;
            }
        }
        maxes[i] = maxpos;
        int columns = wdth - 1;
        int width = (maxpos%mod)*columns/mod;
        for (int j = 0; j < width; j++) fprintf(stderr," ");
        fprintf(stderr,"%s\n",i%2==0?"\e[31mO\e[0m":"\e[32mX\e[0m");
    }

    free(buffer);

    snd_pcm_close (capture_handle);
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    int n = time*rate/buffer_frames;
    for (i = 0; i < n; ++i)
    {
        y+=maxes[i];
        xx+=i*i;
        x+=i;
        xy+=i*maxes[i];
        yy+=maxes[i]*maxes[i];
    }

    x=x*buffer_frames/rate;
    y=y/rate;
    yy=yy/rate/rate;
    xx=xx*buffer_frames/rate*buffer_frames/rate;
    xy=xy*buffer_frames/rate/rate;
    double a = (y*xx-x*xy)/(n*xx-x*x);
    double b=(n*xy-x*y)/(n*xx-x*x);
    double s = sqrt(( yy -2*a*y-2*b*xy+2*a*b*x+a*a*n+b*b*xx)/n);
    fprintf(stderr,"raw rate: %f s/d\n",-b*86400);
    x = 0.;
    y = 0.;
    xx = 0.;
    xy = 0.;
    yy = 0.;
    int m = 0;
    if (verbose)
    {
        for (i = 0; i < time * rate/buffer_frames; ++i)
        {
            printf("%f %f %f %f\n",
                    (double)i*buffer_frames/rate,(double)maxes[i]/rate,
                    ((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b)),s);
        }
    }

    double e;

    for (i = 0; i < n; ++i)
    {
        e = (((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
        if (e*e < threshold*threshold)
        {
            y+=maxes[i];
            xx+=i*i;
            x+=i;
            xy+=i*maxes[i];
            yy+=maxes[i]*maxes[i];
            m++;
        }
    }
    x=x*buffer_frames/rate;
    y=y/rate;
    yy=yy/rate/rate;
    xx=xx*buffer_frames/rate*buffer_frames/rate;
    xy=xy*buffer_frames/rate/rate;
    a = (y*xx-x*xy)/(m*xx-x*x);
    b=(m*xy-x*y)/(m*xx-x*x);
    s = sqrt(( yy -2*a*y-2*b*xy+2*a*b*x+a*a*m+b*b*xx)/m);
    fprintf(stderr,"after %.1fσ removal: %.2f s/d\n",threshold,-b*86400);
  /*  if (verbose)
    {
        for (i = 0; i < time * rate/buffer_frames; ++i)
        {
            e = (((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
            if (e*e < threshold*threshold)
            {
            printf("%f %f %f\n",
                    (double)i*buffer_frames/rate,(double)maxes[i]/rate,
                    ((double)maxes[i]/rate-( a+(double)i*buffer_frames/rate*b))/s);
            }
        }
    }
*/
    exit (0);
}
