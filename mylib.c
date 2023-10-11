#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

snd_pcm_t * initAudio(snd_pcm_format_t format, char* device, unsigned int rate)
{
    int err;
    snd_pcm_t *capture_handle;
    snd_pcm_hw_params_t *hw_params;

    if ((err = snd_pcm_open (&capture_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0) {
        fprintf (stderr, "cannot open audio device %s (%s)\n",
                device,
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
    return capture_handle;
}


fftw_complex * makeFilter(int evalue, int buffer_frames)
{
    fftw_complex *in2 = fftw_alloc_complex(buffer_frames);
    fftw_complex *filterFFT = fftw_alloc_complex(buffer_frames);
    fftw_plan makefilter = fftw_plan_dft_1d(buffer_frames, in2,  filterFFT, FFTW_FORWARD,  FFTW_ESTIMATE);

    // make filter array
    in2[0][0] =1.0;
    for (int j=1; j < buffer_frames; j++)
    {
        in2[j][0] = (evalue==0)?0.0:.398942280401/evalue*(exp(-((float)(j*j))/(float)(evalue*evalue)/2) 
                + exp(-((float)(buffer_frames-j)*(buffer_frames-j))/(float)(evalue*evalue)/2));
        in2[j][1] = 0.0;
    }

    fftw_execute(makefilter);
    fftw_destroy_plan(makefilter);
    fftw_free(in2);
    filterFFT[0][0] = 0.0;
    return filterFFT;
}


int fftfit(int *mean, int *total, int *base, int *val, const fftw_complex *filterFFT, int buffer_frames)
{
    fftw_complex *in = fftw_alloc_complex(buffer_frames);
    fftw_complex *in2 = fftw_alloc_complex(buffer_frames);
    fftw_complex *out = fftw_alloc_complex(buffer_frames);
    fftw_complex *conv = fftw_alloc_complex(buffer_frames);
    fftw_complex *tmp = fftw_alloc_complex(buffer_frames);
    fftw_complex *corr = fftw_alloc_complex(buffer_frames);

    fftw_plan forward = fftw_plan_dft_1d(buffer_frames, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d(buffer_frames, conv, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan corforward = fftw_plan_dft_1d(buffer_frames, in2, tmp, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan correverse = fftw_plan_dft_1d(buffer_frames, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int j=0; j < buffer_frames; j++)
    {
        in[j][0]= (float)mean[j];
        in[j][1] = 0.0;
    }

    fftw_execute(forward);

    for (int j=0; j < buffer_frames ; j++)
    {
        conv[j][0] = (out[j][0]*filterFFT[j][0] - out[j][1]*filterFFT[j][1])/buffer_frames;
        conv[j][1] = (out[j][0]*filterFFT[j][1] + out[j][1]*filterFFT[j][0])/buffer_frames;
    }

    fftw_execute(reverse);

    double ix = 0.0;
    double ixx =0.0;
    double i2x = 0.0;
    double i2xx =0.0;
    for (int j=0; j < buffer_frames ; j++)
    {
        in[j][0] = (float)(in[j][0]/buffer_frames);
        in[j][1] = 0.0;
        ix+=in[j][0];
        in2[j][0] = base[j];
        in2[j][1] = 0.0;
        i2x+=in2[j][0];
    }
    float m=ix/buffer_frames;
    float m2=i2x/buffer_frames;

    for (int j=0; j < buffer_frames ; j++)
    {
        in[j][0] = (in[j][0] - m);
        ixx+=in[j][0]*in[j][0];
        in2[j][0] = (in2[j][0] - m2);
        i2xx+=in2[j][0]*in2[j][0];
    }
    double s = sqrt(ixx/buffer_frames-m/buffer_frames*m/buffer_frames);
    double s2 = sqrt(i2xx*buffer_frames-m2*m2)/buffer_frames;
    // into out
    fftw_execute(forward);

    // into tmp
    fftw_execute(corforward);
    // calculate cross correlation
    for (int j=0; j < buffer_frames ; j++)
    {
        float tmpbuf= tmp[j][0];
        tmp[j][0] = (out[j][0]*tmpbuf + out[j][1]*tmp[j][1])/buffer_frames/buffer_frames/s/s2;
        tmp[j][1] = (-out[j][0]*tmp[j][1] + out[j][1]*tmpbuf)/buffer_frames/buffer_frames/s/s2;
    }
    // transform back into corr
    fftw_execute(correverse);

    float maxcor=-1;
    int poscor=0;
    for (int j=0; j < buffer_frames ; j++)
    {
        if (corr[j][0]>maxcor)
        {
            maxcor =corr[j][0];
            poscor=(j+buffer_frames/2)%buffer_frames;
        }
    }
    // for hexadecimal print 
    *val = (int)(maxcor*16);


    // rescale if large
    if (total[buffer_frames/2]>100000000||total[0]>100)
    {

        long int avg = 0;

        for (int j=0; j < buffer_frames ; j++)
        {
            avg += total[j];

        }
        avg /= buffer_frames;
        int avi = (int)avg;
        if (avi > 100)
        {
            for (int j=0; j < buffer_frames ; j++)
            {
                total[j] -= avi;
            }
        }
        else
        {
            for (int j=0; j < buffer_frames ; j++)
            {
                total[j] /= 2;
            }
        }
    }

    // weigh with square of correlation
    for (int j=0; j < buffer_frames ; j++)
    {
        total[j] = (total[j]+(int)(20*maxcor*maxcor) * mean[(j+poscor+buffer_frames/2+buffer_frames)%buffer_frames]);
    }
    fftw_free(*in);
    fftw_free(*in2);
    fftw_free(*out);
    fftw_free(*conv);
    fftw_free(*tmp);
    fftw_free(*corr);

    return poscor;
}


void linreg(const int *xarr, const int *yarr, int NN, double *a, double *b, double *s)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (int i = 0; i < NN; ++i)
    {
        y  += yarr[i];
        xx += xarr[i]*xarr[i];
        x  += xarr[i];
        xy += xarr[i]*yarr[i];
        yy += yarr[i]*yarr[i];
    }
    
    *a = (y*xx-x*xy)/(NN*xx-x*x);
    *b = (NN*xy-x*y)/(NN*xx-x*x);
    *s = sqrt(( yy -2*(*a)*y-2*(*b)*xy+2*(*a)*(*b)*x+(*a)*(*a)*NN+(*b)*(*b)*xx)/NN);
}
