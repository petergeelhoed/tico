/**
 * Read and parse a wave file
 * PG
* assumes 6Hz and 48000 wav
 **/
#include <unistd.h>
#include <math.h>
#include <fftw3.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0
#define DEBUG 0

// WAVE header structure
 // WAVE file header format
struct HEADER {
    unsigned char riff[4];                      // RIFF string
    unsigned int overall_size   ;               // overall size of file in bytes
    unsigned char wave[4];                      // WAVE string
    unsigned char fmt_chunk_marker[4];          // fmt string with trailing null char
    unsigned int length_of_fmt;                 // length of the format data
    unsigned int format_type;                   // format type. 1-PCM, 3- IEEE float, 6 - 8bit A law, 7 - 8bit mu law
    unsigned int channels;                      // no.of channels
    unsigned int sample_rate;                   // sampling rate (blocks per second)
    unsigned int byterate;                      // SampleRate * NumChannels * BitsPerSample/8
    unsigned int block_align;                   // NumChannels * BitsPerSample/8
    unsigned int bits_per_sample;               // bits per sample, 8- 8bits, 16- 16 bits etc
    unsigned char data_chunk_header [4];        // DATA string or FLLR string
    unsigned int data_size;                     // NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
};


unsigned char buffer4[4];
unsigned char buffer2[2];

struct HEADER header;

int main(int argc, char **argv) 
{
    FILE *outfile=stdout;
    int uvalue = 32768;
    int nvalue = 48000;
    int dvalue = nvalue;
    int c;
    int cvalue = 0;
    int evalue = 0;
    int ovalue = 0;
    int pvalue = 0;
    int fvalue = 48000;
    int hvalue = 21600;
    int lvalue = 0;
    int rvalue = 0;
    int qvalue = 2000;
    int mvalue = 0;
    int svalue = nvalue;
    int tvalue = 0;
    int wvalue = 0;
    int vvalue = 0;
    int NN;
    opterr=0;

    while ((c = getopt (argc, argv, "c:u:n:d:l:r:q:m:s:twvh:f:e:op:")) != -1)
        switch (c)
        {
            case 'u':
                // upper limit
                uvalue = atoi(optarg);
                break;
            case 'n':
                // maximum number of peak points
                nvalue = atoi(optarg);
                break;
            case 's':
                // around m
                svalue = atoi(optarg);
                mvalue=0?-1:mvalue;
                break;
            case 'p':
                //print teeth hisdev
                pvalue = atoi(optarg)*2;
                fprintf(stderr,"teeth: %d\n",pvalue/2);
                wvalue = 1;
                ovalue = 1;
                break;
            case 'm':
                //middle
                mvalue = atoi(optarg);
                break;
            case 'o':
                //mean output normalised and shifted
                wvalue = 1;
                ovalue = 1;
                break;
            case 'v':
                //verbose
                vvalue = 1;
                break;
            case 'e':
                //exp filter s=optarg 2stdev
                evalue = atoi(optarg);
                break;
            case 'w':
                //raw output
                wvalue = 1;
                break;
            case 't':
                //peak only
                tvalue = 1;
                break;
            case 'h':
                //trim left
                hvalue = atoi(optarg);
                break;
            case 'f':
                //trim left
                fvalue = atoi(optarg);
                break;
            case 'l':
                //trim left
                lvalue = atoi(optarg);
                break;
            case 'q':
                //move up 
                qvalue = atoi(optarg);
                break;
            case 'r':
                //trim right
                rvalue = atoi(optarg);
                break;
            case 'd':
                //maximum distandce to peak 
                dvalue = atoi(optarg);
                break;
            case 'c':
                // lower limit
                cvalue = atoi(optarg);
                break;
            case '?':
                if (optopt == 'r' || optopt == 'l' || optopt == 'c' || optopt == 'u' || optopt == 'n' || optopt == 'd' || optopt == 'm')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt)){
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    fprintf (stderr, "h bph default 21600\nf wav frequency default 48000\nc lower limit\nu upper limit\nn maximum points\nd max distance\nl left trim (s)\nq move points up (default 2000)\nr right trim (s)\nm mean, use with -s\nt show only peakvalue\nw raw input\nv show gnuplot command\ne gausfiliter stdev\np teethfor hisdev\n");
                }else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort ();
        }

    NN=(int)(fvalue*3600/hvalue);
    lvalue = lvalue*hvalue/3600;
    rvalue = rvalue*hvalue/3600;
    if (DEBUG != 0) fprintf(stderr,"Cutoff  %d < value < %d  ; n < %d ; d < %d; m = %d+-%d\n",cvalue, uvalue,nvalue,dvalue ,mvalue,svalue);


    outfile = fopen("oink", "w");


    int read = 0;

    // read header parts

    read = fread(header.riff, sizeof(header.riff), 1, stdin);
    if (DEBUG != 0) printf("(1-4): %s \n", header.riff);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0) printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    // convert little endian to big endian 4 byte int
    header.overall_size  = buffer4[0] |
        (buffer4[1]<<8) |
        (buffer4[2]<<16) |
        (buffer4[3]<<24);

    if (DEBUG != 0) printf("(5-8) Overall size: bytes:%u, Kb:%u \n", header.overall_size, header.overall_size/1024);

    read = fread(header.wave, sizeof(header.wave), 1, stdin);
    if (DEBUG != 0) printf("(9-12) Wave marker: %s\n", header.wave);

    read = fread(header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, stdin);
    if (DEBUG != 0) printf("(13-16) Fmt marker: %s\n", header.fmt_chunk_marker);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0) printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    // convert little endian to big endian 4 byte integer
    header.length_of_fmt = buffer4[0] |
        (buffer4[1] << 8) |
        (buffer4[2] << 16) |
        (buffer4[3] << 24);
    if (DEBUG != 0) printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);

    read = fread(buffer2, sizeof(buffer2), 1, stdin); if (DEBUG != 0) printf("%u %u \n", buffer2[0], buffer2[1]);

    header.format_type = buffer2[0] | (buffer2[1] << 8);
    char format_name[10] = "";
    if (header.format_type == 1)
        strcpy(format_name,"PCM");
    else if (header.format_type == 6)
        strcpy(format_name, "A-law");
	 else if (header.format_type == 7)
	  strcpy(format_name, "Mu-law");
	
	 if (DEBUG != 0) printf("(21-22) Format type: %u %s \n", header.format_type, format_name);
	
	 read = fread(buffer2, sizeof(buffer2), 1, stdin);
	 if (DEBUG != 0) printf("%u %u \n", buffer2[0], buffer2[1]);
	
	 header.channels = buffer2[0] | (buffer2[1] << 8);
	 if (DEBUG != 0) printf("(23-24) Channels: %u \n", header.channels);
	
	 read = fread(buffer4, sizeof(buffer4), 1, stdin);
	 if (DEBUG != 0) printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
	
	 header.sample_rate = buffer4[0] |
	                        (buffer4[1] << 8) |
	                        (buffer4[2] << 16) |
	                        (buffer4[3] << 24);
	
	 if (DEBUG != 0) printf("(25-28) Sample rate: %u\n", header.sample_rate);
	
	 read = fread(buffer4, sizeof(buffer4), 1, stdin);
	 if (DEBUG != 0) printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
	
	 header.byterate  = buffer4[0] |
	                        (buffer4[1] << 8) |
	                        (buffer4[2] << 16) |
	                        (buffer4[3] << 24);
	 if (DEBUG != 0) printf("(29-32) Byte Rate: %u , Bit Rate:%u\n", header.byterate, header.byterate*8);
	
	 read = fread(buffer2, sizeof(buffer2), 1, stdin);
	 if (DEBUG != 0) printf("%u %u \n", buffer2[0], buffer2[1]);
	
	 header.block_align = buffer2[0] |
	                    (buffer2[1] << 8);
	 if (DEBUG != 0) printf("(33-34) Block Alignment: %u \n", header.block_align);
	
	 read = fread(buffer2, sizeof(buffer2), 1, stdin);
	 if (DEBUG != 0) printf("%u %u \n", buffer2[0], buffer2[1]);
	
	 header.bits_per_sample = buffer2[0] |
	                    (buffer2[1] << 8);
	 if (DEBUG != 0) printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);
	
	 read = fread(header.data_chunk_header, sizeof(header.data_chunk_header), 1, stdin);
	 if (DEBUG != 0) printf("(37-40) Data Marker: %s \n", header.data_chunk_header);
	
	 read = fread(buffer4, sizeof(buffer4), 1, stdin);
	 if (DEBUG != 0) printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
	
	 header.data_size = buffer4[0] |
	                (buffer4[1] << 8) |
	                (buffer4[2] << 16) |
	                (buffer4[3] << 24 );
	 if (DEBUG != 0) printf("(41-44) Size of data chunk: %u \n", header.data_size);
	
	
	 // calculate no.of samples
	 long num_samples = (8 * header.data_size) / (header.channels * header.bits_per_sample);
	 if (DEBUG != 0) printf("Number of samples:%lu \n", num_samples);
	
	 long size_of_each_sample = (header.channels * header.bits_per_sample) / 8;
	 if (DEBUG != 0) printf("Size of each sample:%ld bytes\n", size_of_each_sample);
	
	 // calculate duration of file
	 float duration_in_seconds = (float) header.overall_size / header.byterate;
	 if (DEBUG != 0) printf("Approx.Duration in seconds=%f\n", duration_in_seconds);
	
	
	
	 // read each sample from data chunk if PCM
	 if (header.format_type == 1) { // PCM
         long i =0;
         unsigned char lsb[1];
         signed char msb[1];

         for (i =44; i <= 100 ; i++) {
             //skip shit
             read = fread(lsb, sizeof(lsb), 1, stdin);
             read += fread(msb, sizeof(msb), 1, stdin);
         }
         int mean[NN];
         int peaks[NN];

         float x[NN*2][pvalue];
         float xx[NN*2][pvalue];
         int nn[NN*2][pvalue];

         fftw_complex *in, *out, *out2, *conv,  *in2; /* double [2] */
         fftw_plan p, q, pr;

         in = fftw_alloc_complex(NN);
         in2 = fftw_alloc_complex(NN);
         for (int j=0; j < NN; j++) 
         {
             in2[j][0] = .398942280401/evalue*(exp(-((float)(j*j))/(float)(evalue*evalue)/2) + exp(-((float)(NN-j)*(NN-j))/(float)(evalue*evalue)/2));
             //in2[j][0] = (exp(-((float)(j*j))/(float)(evalue*evalue)/2) + exp(-((float)(NN-j)*(NN-j))/(float)(evalue*evalue)/2));
             in2[j][1] = 0.0;
         }
         out = fftw_alloc_complex(NN);
         out2 = fftw_alloc_complex(NN);
         conv = fftw_alloc_complex(NN);

         // filter
         p = fftw_plan_dft_1d(NN, in,out, FFTW_FORWARD, FFTW_ESTIMATE );
         q = fftw_plan_dft_1d(NN, in2,out2, FFTW_FORWARD, FFTW_ESTIMATE );
         pr = fftw_plan_dft_1d(NN, conv, in, FFTW_BACKWARD, FFTW_ESTIMATE);
         fftw_execute(q);
         fftw_destroy_plan(q);

         int val = 0;
         for (i =qvalue+1; i <= header.overall_size/2+qvalue-rvalue*NN-100 ; i++) 
         {
             read = fread(lsb, sizeof(lsb), 1, stdin);
             read += fread(msb, sizeof(msb), 1, stdin);
             if (read == 2) {
                 val = 0 ;
                 mean[i%NN] = 0;
                 val = abs((msb[0] << 8) | lsb[0]);
                 if (val >= cvalue && val <= uvalue ){
                     mean[i%NN] = val ; 
                 }
                 //at the end of the loop:
                 if (i%NN == 0)
                 {
                     // filter mean 
                     if (evalue > 0)
                     {
                         for (int j=0; j < NN; j++) 
                         {
                             in[j][0] = (float)mean[j];
                             in[j][1] = 0.0;
                         }
                         fftw_execute(p);

                         for (int j=0; j < NN ; j++)
                         {
                             conv[j][0] = (
                                     +out[j][0]*out2[j][0]
                                     -out[j][1]*out2[j][1]);
                             conv[j][1] = (
                                     out[j][0]*out2[j][1]
                                     +out[j][1]*out2[j][0]);

                         }


                         fftw_execute(pr);


                         for (int j=0; j < NN ; j++)
                         {
                             mean[j] = (int)(in[j][0]/NN);
                         }

                         
                     }
                     
                     float tot=0;
                     float mom=0;
                     int center = 0;
                     int j=0;
                     int n=0;
                     int max=-1;
                     int k=0;
                     for (j=0; j < NN ; j++) 
                     {
                         n+=(mean[j]>0);
                         tot += mean[j];
                         mom += j*mean[j];
                         if (mean[j]>max && mean[j]>cvalue)
                         {
                             max = mean[j];
                             k = j;
                         }
                     }
                     peaks[(int)(i-1)/NN]=(max>0&&abs(k-mvalue < svalue ))?k:-1;
                     if (wvalue && i>=lvalue*NN)
                     {
                         for (int j=0; j < NN ; j++)
                         {
                             if (ovalue == 0)
                             {
                                 printf("%d %d %d\n",i+j-NN,mean[j],peaks[(int)(i-1)/NN]);
                             }else{
                                 if (pvalue == 0) 
                                 {
                                     printf("%d %f %d\n",j-peaks[(int)(i-1)/NN],(float)mean[j]/max,(int)((i-1)/NN));
                                 }else{ 
                                     x[j-peaks[(int)(i-1)/NN]+NN][(int)(i-1)/NN%pvalue] 
                                         += (float)mean[j]/max;
                                     xx[j-peaks[(int)(i-1)/NN]+NN][(int)(i-1)/NN%pvalue] 
                                         += (float)mean[j]*mean[j]/max/max;
                                     nn[j-peaks[(int)(i-1)/NN]+NN][(int)(i-1)/NN%pvalue] += 1;
                                 }

                             }
                         }
                     }

                     if (tvalue > 0 ) 
                     {
                         if((i>=lvalue*NN && n <= nvalue && max > 0 &&
                                     ((mvalue > 0  && abs(k-mvalue) < svalue) || mvalue == 0)))  {
                             fprintf(outfile,"%d %d %d\n",(i-1)/NN,k,max);
                         } 
                     }
                     else if (tot>0) {
                         center = (int)(mom/tot);
                         for (j=0; j < NN ; j++) 
                         {
                             if (i>=lvalue*NN && mean[j] > 0 && n <= nvalue &&
                                     (abs(j-center) < dvalue) &&
                                     ((mvalue > 0  && abs(j-mvalue) < svalue) || mvalue == 0))  {
                                 fprintf(outfile,"%d %d %d\n",(i-1)/NN,j,mean[j]);

                             }
                         }
                     }

                 }
             } else {
                 printf("Error reading file. %d bytes\n", read);
                 exit -1;
             }
         }
         fftw_destroy_plan(p);
         fftw_destroy_plan(pr);
         fftw_cleanup();
             for (int k=0; k < pvalue ; k++)
             {
             for (int j=0; j < 2*NN; j++)
             {
                 if (nn[j][k] > 1 )
                 {
                     printf("%d %f %f %d %d\n",
                             j-NN,
                             x[j][k]/nn[j][k]*(1-2*(k%2==0)),
                             sqrt(xx[j][k]/nn[j][k]-x[j][k]*x[j][k]/nn[j][k]/nn[j][k]),nn[j][k],
                             (int)(((int)(k+2)/2)*(1-2*(k%2==0))));
                 }
             }
             printf("\n");
             }
         float y=0;
         float yy=0;
         int n=0;
         for (int j=0; j < (i-1)/NN ; j++) 
         {
             if (peaks[j]>0)
             {
                 n++;
                 y+=peaks[j];
                 yy+=peaks[j]*peaks[j];
             }
         }
         fprintf(stderr,"points: %d  mean:%f  stdev:%f\n",n,y/n,sqrt(yy/n-y/n*y/n));
     }
     fclose(outfile);

     char command[1024] ;
     sprintf(command," echo 'uns colorbox; f=%d;h=%d/3600; set cbtics 1; set term png size 1920,1080 font \"DejaVuSansCondensed,12 truecolor \" ; set ytics nomirror; set out \"/dev/null\"; plot \"oink\" u ($1/h):($2/f); set y2tics ; set y2range [GPVAL_Y_MIN*f:GPVAL_Y_MAX*f]; set out \"/home/pi/lussen/www/tico.png\"; set fit quiet; set samples 1000; d=15;d1=15; f(x)=(a-b*x+c*cos(3.1415926*(x-x0)*h/d)); g(x)=(a1-b1*x+c1*cos(3.1415926*(x-x1)*h/d1)); fit []f(x) \"oink\" u ($1/h):(int($1)%%2==0?$2/f:NaN) via a,b,c,d,x0; print d,c*1000,\"ms \",b*86400,\"s/d\"; fit []g(x) \"oink\" u ($1/h):(int($1)%%2==1?$2/f:NaN) via a1,b1,c1,d1,x1; print d1,c1*1000,\"ms \",b1*86400,\"s/d\"; set xrange [:]; set key below;set samples 1000; set xlabel \"time (s)\"; set ylabel sprintf(\"modulo 1/%d s (s)\",h);plot \"oink\" u ($1/h):($2/f):(int($1)%%2) pal ps 3 pt  5 t sprintf(\"%%.1f s/d %%.1f s/d\",b*86400,b1*86400) , f(x) lc 7 t sprintf(\"%%.2fms\",c*1000), g(x) lc 5 t sprintf(\"%%.2fms\",c1*1000) ; print \"beaterror: \",1000*(a-a1), \"ms\";'  | gnuplot -persist ",fvalue,hvalue,hvalue/3600);


     if (vvalue) fprintf(stderr,"%s\n",command);
     return system(command);
}
