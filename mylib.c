#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>


void normalise(int NN,fftw_complex *in)
{
    double ix = 0.0;
    double ixx =0.0;
    for (int j = 0; j < NN ; j++)
    {
        in[j][0] = in[j][0];
        ix+=in[j][0];
        ixx+=in[j][0]*in[j][0];
    }
    double m = ix/NN;
    double s = sqrt(ixx/NN-m*m);
    for (int j = 0; j < NN ; j++)
    {
        in[j][0] = (in[j][0]-m)/s;
    }
}


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


fftw_complex * makeFilter(int evalue, int NN)
{
    fftw_complex *in2 = fftw_alloc_complex(NN);
    fftw_complex *filterFFT = fftw_alloc_complex(NN);
    fftw_plan makefilter = fftw_plan_dft_1d(NN, in2,  filterFFT, FFTW_FORWARD,  FFTW_ESTIMATE);

    if (evalue != 0)
    {
        // make filter array
        for (int j = 0; j < evalue*5; j++)
        {
            in2[j][0] = .398942280401/evalue*( exp(-((double)(j*j))/(double)(evalue*evalue)/2) );
            in2[j][1] = 0.0;
        }
        for (int j = evalue*5 ; j < NN-evalue*5; j++)
        {
            in2[j][0] = 0.0;
            in2[j][1] = 0.0;
        }

        for (int j = NN- evalue*5 ; j < NN; j++)
        {
            in2[j][0] = .398942280401/evalue*( exp(-((double)(NN-j)*(NN-j))/(double)(evalue*evalue)/2));
            in2[j][1] = 0.0;
        }
    }
    else
    {
        in2[0][0] = 100; 
        in2[0][1] = 0; 
        for (int j = 1; j < NN; j++)
        {
            in2[j][0] = 0; 
            in2[j][1] = 0; 
        }
    }
    

    fftw_execute(makefilter);
    fftw_destroy_plan(makefilter);
    fftw_free(in2);
    return filterFFT;
}


fftw_complex* convolute(int NN, int* array, const fftw_complex *filterFFT)
{
    fftw_complex *in = fftw_alloc_complex(NN);
    fftw_complex *out = fftw_alloc_complex(NN);
    fftw_plan forward = fftw_plan_dft_1d(NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d(NN, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int j = 0; j < NN; j++)
    {
        in[j][0] = (double)array[j];
        in[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (int j = 0; j < NN ; j++)
    {
        out[j][0] = (out[j][0]*filterFFT[j][0] - out[j][1]*filterFFT[j][1])/NN;
        out[j][1] = (out[j][0]*filterFFT[j][1] + out[j][1]*filterFFT[j][0])/NN;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);
    fftw_free(*out);
    return in;
}


void rescale(int* total, int NN)
{
    if (total[NN/2]>100000000||total[0]>100)
    {

        long int avg = 0;

        for (int j = 0; j < NN ; j++)
        {
            avg += total[j];

        }
        avg /= NN;
        int avi = (int)avg;
        if (avi > 100)
        {
            for (int j = 0; j < NN ; j++)
            {
                total[j] -= avi;
            }
        }
        else
        {
            for (int j = 0; j < NN ; j++)
            {
                total[j] /= 2;
            }
        }
    }
}


fftw_complex* crosscor(int NN, fftw_complex* array, fftw_complex* ref) 
{
    normalise(NN, array);
    normalise(NN, ref);
    fftw_complex *tmparr = fftw_alloc_complex(NN);
    fftw_complex *tmpref = fftw_alloc_complex(NN);
    fftw_complex *tmp = fftw_alloc_complex(NN);
    fftw_complex *corr = fftw_alloc_complex(NN);

    fftw_plan arrFour = fftw_plan_dft_1d(NN, array, tmparr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan refFour= fftw_plan_dft_1d(NN, ref, tmpref, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan correverse = fftw_plan_dft_1d(NN, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(arrFour);
    fftw_execute(refFour);

    for (int j = 0; j < NN ; j++)
    {
        tmp[j][0] = (tmparr[j][0]*tmpref[j][0] + tmparr[j][1]*tmpref[j][1]);
        tmp[j][1] = (-tmparr[j][0]*tmpref[j][1] + tmparr[j][1]*tmpref[j][0]);
    }
    fftw_execute(correverse);

    double scale = 1./NN/NN;
    for (int j = 0; j < NN ; j++)
    {
        corr[j][0] *= scale;
    }

    fftw_destroy_plan(correverse);
    fftw_destroy_plan(arrFour);
    fftw_destroy_plan(refFour);

    fftw_free(*tmparr);
    fftw_free(*tmpref);
    fftw_free(*tmp);
    return corr;
}


void applyFilter(int* input, int NN, fftw_complex* filterFFT, double* out)
{
    fftw_complex *filteredinput = convolute(NN,input,filterFFT);
    for (int j = 0; j < NN ; j++)
    {
        out[j]=filteredinput[j][0];
    }
    fftw_free(filteredinput);
}


int fftfit(
        int* input,
        int* total,
        int* base,
        int* hexvalue,
        const fftw_complex *filterFFT,
        int NN,
        int halfsearch)
{
    fftw_complex *Fbase = fftw_alloc_complex(NN);
    fftw_complex *filteredinput = convolute(NN,input,filterFFT);


    for (int j = 0; j < NN ; j++)
    {
        Fbase[j][0] = (double)base[j];
        Fbase[j][1] = 0.0;
    }

    fftw_complex* corr = crosscor(NN,filteredinput,Fbase);

    double maxcor = -1;
    int poscor = 0;
    int startsearch = halfsearch?(NN/4):0;
    int stopsearch = halfsearch?(NN*3/4):NN;
    for (int j = startsearch; j < stopsearch ; j++)
    {
        if (corr[j][0]>maxcor)
        {
            maxcor =corr[j][0];
            poscor = j;
        }
    }

    // for hexadecimal print 
    *hexvalue = (int)(maxcor*16);

    poscor -= NN/2;

    // rescale if large
    if (total)
    {
        rescale(total, NN);

        // weigh with square of correlation
        for (int j = 0; j < NN ; j++)
        {
            total[j] = (total[j]+(int)(2000*maxcor*maxcor) * filteredinput[(j+poscor+NN/2+NN)%NN][0]);
        }
    }
    fftw_free(*filteredinput);
    fftw_free(*Fbase);
    fftw_free(*corr);

    return poscor;
}


void readBufferRaw(snd_pcm_t *capture_handle, int NN, char *buffer, int* in)
{
        unsigned char lsb;
        signed char msb;
        int err;
        if ((err = snd_pcm_readi (capture_handle, buffer, NN)) != NN) 
        {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
        for (int j = 0; j < NN*2; j+=2) 
        {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2] = (msb << 8) | lsb ;
        }
}


void readBuffer(snd_pcm_t *capture_handle, int NN, char *buffer, int* derivative)
{
        int in[NN];
        unsigned char lsb;
        signed char msb;
        int err;
        if ((err = snd_pcm_readi (capture_handle, buffer, NN)) != NN) 
        {
            fprintf (stderr, "read from audio interface failed %d (%s)\n", err, snd_strerror (err));
            exit (1);
        }
        for (int j = 0; j < NN*2; j+=2) 
        {
            msb = *(buffer+j+1);
            lsb = *(buffer+j);
            in[j/2] = (msb << 8) | lsb ;

            derivative[j/2] = (j==0)?0:fabs(in[j/2]-in[j/2-1]);
        }
}


void readShiftedBuffer(int* derivative, snd_pcm_t *capture_handle, int NN, char* buffer, int maxpos, int shift, int* totalshift, int lowerBound, int upperBound)
{

    if (maxpos < lowerBound)
    {
        *totalshift -= shift;
        memcpy(derivative+NN-shift, derivative , shift*sizeof(int));
        readBuffer(capture_handle, NN-shift, buffer, derivative);
    }
    else if (maxpos > upperBound ) 
    {
        *totalshift += shift;
        readBuffer(capture_handle, shift, buffer, derivative);
        readBuffer(capture_handle, NN, buffer, derivative);
    }
    else
    {
        readBuffer(capture_handle, NN, buffer, derivative);
    }
}


void printspaces(int maxpos,int hexvalue, char* spaces,int mod,int columns, double a,double b,int NN,int i,float beatError)
{
    while (maxpos < mod) maxpos+=mod;
    while (a < mod) a+=(double)mod;

    int width = (maxpos%mod)*columns/mod;
    int widtha = (((int)a)%mod)*columns/mod;
    fprintf(stderr,"%4.1fms %6.1fs/d",beatError,b*86400/NN);
    memset(spaces, ' ', columns);
    spaces[widtha] = '|';
    spaces[width] = '\0';
    fprintf(stderr,"%s%s%X\e[0m",spaces,i%2==0?"\e[31m": "\e[32m",hexvalue);
    memset(spaces, ' ', columns);
    if (widtha > width)
    {
        spaces[widtha-width-1] = '|';
        spaces[widtha-width-1+1] = '\0';
        fprintf(stderr,"%s",spaces);
    }
    fprintf(stderr,"\n");
}


void linreg(const int* xarr, const int* yarr, int NN, double* a, double* b, double* s)
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
    *s = sqrt((yy-2*(*a)*y-2*(*b)*xy+2*(*a)*(*b)*x+(*a)*(*a)*NN+(*b)*(*b)*xx)/NN);
}


void fit10secs(
        double* a,
        double* b,
        double* s,
        int i,
        int* maxvals,
        int* maxes,
        int qvalue,
        int cvalue,
        int npeaks)
{
    int m = 0;
    int fitwindow = i>npeaks*(1+qvalue)?npeaks*(1+qvalue):i;

    if (i >= fitwindow)
    {
        int xarr[fitwindow];
        int yarr[fitwindow];
        for (int k = 0; k < fitwindow;k+=qvalue+1)
        {
            if (maxvals[i-k] > cvalue)
            {
                yarr[m] = maxes[i-k];
                xarr[m] = k;
                m++;
            }
        }
        if (m > 1)
        {
            linreg(xarr,yarr, m, a, b, s);
        }
    }
}


void writefiles(
        FILE* fptotal,
        FILE* rawfile,
        int* totaltick,
        int* totaltock,
        int* defaultpulse,
        int* maxpos,
        int n,
        int NN)
{
    if (fptotal)
    {
        if (NN == 8000 || NN== 16000)
        {
            for (int j = 0; j < NN; j++) fprintf(fptotal,"%d %d %d\n",totaltick[j],totaltock[j],defaultpulse[j]);
        }
        else
        {
            for (int j = 0; j < NN; j++) fprintf(fptotal,"%d %d\n",totaltick[j],totaltock[j]);
        }
        fclose(fptotal);
    }
    if (rawfile)
    {
        for (int i = 0; i < n ; ++i) fprintf(rawfile,"%d %d\n",i,maxpos[i]);
        fclose(rawfile);
    }
}


void calculateTotal(int n, int* maxpos,int NN, double threshold)
{
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int xarr[n];

    for (int i = 0; i < n ; ++i)
    {
        xarr[i] = i;
    }

    linreg(xarr,maxpos, n, &a, &b, &s);

    /*
       a /= NN*NN;
       b /= NN;
       s /= rate;
     */

    fprintf(stderr,"raw rate: %f s/d\n",-b*86400/NN);
    int m = 0;

    double e;

    for (int i = 0; i < n; ++i)
    {
        e = fabs(((double)maxpos[i]-(a+xarr[i]*b))/s);
        if (e < threshold)
        {
            maxpos[m] = maxpos[i];
            xarr[m] = xarr[i];
            m++;
        }
    }
    linreg(xarr, maxpos, m, &a, &b, &s);

    fprintf(stderr,"after %.1fσ removal: %.2f s/d\n",threshold,-b*86400/NN);
}

int getBeatError(int* totaltick, int NN)
{
    int maxtick = -1;
    int postick = 0;
    for (int j = 0; j < NN/2; j++)
    {
        if ( totaltick[j] > maxtick )
        {
            maxtick = totaltick[j];
            postick = j; 
        }
    }
    int maxtock = -1;
    int postock = NN/2;
    for (int j = NN/2 ; j< NN ; j++)
    {
        if ( totaltick[j] > maxtock )
        {
            maxtock = totaltick[j];
            postock = j; 
        }
    }
    return postock-postick-NN/2;

}

void crosscorint(int NN, int* array, int* ref, int* cross ) 
{

    fftw_complex *tmparr = fftw_alloc_complex(NN);
    fftw_complex *tmpref = fftw_alloc_complex(NN);
    for (int j = 0; j < NN ; j++)
    {
        tmparr[j][0] = array[j];
        tmparr[j][1] = 0.0;
        tmpref[j][0] = ref[j];
        tmpref[j][1] = 0.0;
    }
    fftw_complex* coor =  crosscor(NN, tmparr, tmpref);

    for (int j = 0; j < NN ; j++)
    {
        cross[j] = (int)(coor[j][0]*NN);
    }

}
