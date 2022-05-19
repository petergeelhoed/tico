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

struct HEADER readheader();
int openfiles(FILE **tickfile,  FILE **tockfile, 
              FILE **corfile,   FILE **rawfile, 
              FILE *pulsefile, FILE *pulsefiletock,
			  double *ps,      double *pst,
			  int jvalue, int wvalue, int NN);

int main(int argc, char **argv) 
{
    FILE *tickfile=stdout;
    FILE *tockfile=stdout;
    FILE *corfile=stdout;
    FILE *rawfile=stdout;
    FILE *pulsefile;
    FILE *pulsefiletock;
    int nvalue = 48000;
    int dvalue = 4000;
    int c;
    int evalue = 4;
    float cvalue = 0;
    int ovalue = 0;
    int pvalue = 15*2;
    int fvalue = 48000;
    int hvalue = 21600;
    int lvalue = 0;
    int rvalue = 0;
    int qvalue = 4000;
    int tvalue = 0;
    int kvalue = 0;
    int wvalue = 0;
    int vvalue = 0;
    int jvalue = 0;
    int NN     = 8000;
    opterr=0;
    double *ps;
    double *pst;
    int read = 0;

    while ((c = getopt (argc, argv, "n:d:l:r:q:twvh:f:e:op:jkc:")) != -1)
        switch (c)
        {
            case 'n':
                // maximum number of peak points
                nvalue = atoi(optarg);
                break;
            case 'p':
                //print teeth hisdev
                pvalue = atoi(optarg)*2;
                fprintf(stderr,"teeth: %d\n",pvalue/2);
                break;
            case 'k':
                //no corrshift
                kvalue = 1;
                break;
            case 'j':
                //indata
                jvalue = 1;
                break;
            case 'o':
                //use loudest noise instead of correlationpeak
                ovalue = 1;
                break;
            case 'v':
                //verbose
                vvalue = 1;
                break;
            case 'c':
                //exp filter s=optarg 2stdev
                cvalue = atof(optarg);
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
                // toglle tic tock?
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
                //max windowshift
                dvalue = atoi(optarg);
                break;
            case '?':
                if (optopt == 'r' || optopt == 'l' || optopt == 'c' || optopt == 'u' || optopt == 'n' || optopt == 'd' )
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt)){
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    fprintf (stderr, "h bph default 21600\nf wav frequency default 48000\nn maximum points\nd max distance for window shift\nl left trim (s)\nq move points up (default 2000)\nr right trim (s)\nm mean, use with -s\nj flatten the curve \nw raw input\nv show gnuplot command\ne gausfiliter stdev\np teethfor hisdev\nt toggle tick/tock\n");
                }else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort ();
        }

    NN=(int)(fvalue*3600/hvalue);
    pst = malloc(NN*sizeof(double));
    ps = malloc(NN*sizeof(double));

    lvalue = lvalue*hvalue/3600;
    rvalue = rvalue*hvalue/3600;
	if (openfiles(&tickfile,  &tockfile, 
				  &corfile,   &rawfile, 
				   pulsefile, pulsefiletock,
				   ps, pst, jvalue, wvalue,NN))
	{
		fprintf (stderr,"files errored out\n");
		exit(1);
	}

	struct HEADER header = readheader();
	// read each sample from data chunk if PCM
	if (header.format_type == 1) 
	{ // PCM
         long i =0;
         unsigned char lsb[1];
         signed char msb[1];
		 int globalshift = 0;
         float mean[NN];
		 mean[0]=1.;
		 for(int j=1;j<NN;j++) { mean[j]=0.; }

         for (i =44; i <= 100 ; i++) {
             //skip shit
             read = fread(lsb, sizeof(lsb), 1, stdin);
             read += fread(msb, sizeof(msb), 1, stdin);
         }


         fftw_complex *in, *out, *filterFFT, *conv,  *in2, *tmp,*corr; 
         fftw_plan p, q, pr,cf,cr;

         in = fftw_alloc_complex(NN);
         in2 = fftw_alloc_complex(NN);
         out = fftw_alloc_complex(NN);
         filterFFT = fftw_alloc_complex(NN);
         conv = fftw_alloc_complex(NN);
         tmp = fftw_alloc_complex(NN);
         corr = fftw_alloc_complex(NN);


		 // make filter array
         for (int j=0; j < NN; j++) 
         {
             in2[j][0] = .398942280401/evalue*(exp(-((float)(j*j))/(float)(evalue*evalue)/2) + exp(-((float)(NN-j)*(NN-j))/(float)(evalue*evalue)/2));
             in2[j][1] = 0.0;
         }
         
         p = fftw_plan_dft_1d(NN, in,out, FFTW_FORWARD, FFTW_ESTIMATE );
         q = fftw_plan_dft_1d(NN, in2,filterFFT, FFTW_FORWARD, FFTW_ESTIMATE );
         pr = fftw_plan_dft_1d(NN, conv, in, FFTW_BACKWARD, FFTW_ESTIMATE);
         cf = fftw_plan_dft_1d(NN,  in2,tmp, FFTW_FORWARD, FFTW_ESTIMATE );
         cr = fftw_plan_dft_1d(NN, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);
         fftw_execute(q);
         fftw_destroy_plan(q);

         int val = 0;
         int Npeak =0;
		 int shift = NN;
		 int startshift = NN;
         long length = header.overall_size/2+qvalue-rvalue*NN-100 - NN*2 ; 

		 // skip q values
		 for(int j=0;j<qvalue;j++)
		 {
			 read = fread(lsb, sizeof(lsb), 1, stdin);
			 read += fread(msb, sizeof(msb), 1, stdin);
			 if (read == 2) 
			 {
				 i++;
			 } else {
				 printf("Error reading file. %d bytes\n", read);
				 return -1;
			 }
		 }

		 // loop the entire file
		 while ( i < length) 
		 {
			 // 'reread' the data for smaller shitfs
			 for(int j=shift;j<NN;j++) { mean[j-shift]=mean[j]; }

			 // now get the data you need
             for(int j=NN-shift;j<NN;j++)
			 {
				 read = fread(lsb, sizeof(lsb), 1, stdin);
				 read += fread(msb, sizeof(msb), 1, stdin);
				 if (read == 2) 
				 {
					 i++;
					 mean[j] = (float)abs((msb[0] << 8) | lsb[0]);
				 } else {
					 printf("Error reading file. %d bytes\n", read);
					 exit(-1);
				 }
			 }


			 // filter mean 
			 if (evalue > 0)
			 {
				 for (int j=0; j < NN; j++) 
				 {
					 in[j][0] = mean[j];
					 in[j][1] = 0.0;
				 }
				 fftw_execute(p);

				 for (int j=0; j < NN ; j++)
				 {
					 conv[j][0] = (
							 +out[j][0]*filterFFT[j][0]
							 -out[j][1]*filterFFT[j][1])/NN;
					 conv[j][1] = (
							 out[j][0]*filterFFT[j][1]
							 +out[j][1]*filterFFT[j][0])/NN;
				 }

				 fftw_execute(pr);
				 for (int j=0; j < NN ; j++) { mean[j] = in[j][0]; }
			 }
                     
			 float tot=0;
			 float mom=0;
			 int j=0;
			 int n=0;
			 if (Npeak>= lvalue)
			 {
				 for (int j=0; j < NN ; j++)
				 {
					 in[j][0] = mean[j];
					 in[j][1]= 0.0;
					 if (wvalue) fprintf(rawfile, "%f\n",mean[j]);
				 }

				 double ix = 0.0;
				 double ixx =0.0;
				 for (int j=0; j < NN ; j++)
				 {
					 in[j][0] = (float)(in[j][0]/NN);
					 in[j][1] = 0.0;
					 ix+=in[j][0];
				 }
				 //use in2 for second into tmp
				 double i2x = 0.0;
				 double i2xx =0.0;
				 for (int j=0; j < NN ; j++)
				 {
					 in2[j][0]= (Npeak%2==tvalue)?pst[j]:ps[j];
					 in2[j][1] = 0.0;
					 i2x+=in2[j][0];
				 }
				 float m=ix/NN;
				 float m2=i2x/NN;

				 for (int j=0; j < NN ; j++)
				 {
					 in[j][0] = (in[j][0] - m);
					 ixx+=in[j][0]*in[j][0];
					 in2[j][0] = (in2[j][0] - m2);
					 i2xx+=in2[j][0]*in2[j][0];
				 }
				 double s = sqrt(ixx/NN-m/NN*m/NN);
				 double s2 = sqrt(i2xx*NN-m2*m2)/NN;
				 fftw_execute(p);
				 fftw_execute(cf);
				 // calculate cross correlation
				 for (int j=0; j < NN ; j++)
				 {
					 float tmpbuf= tmp[j][0];
					 tmp[j][0] = (
							 +out[j][0]*tmpbuf
							 +out[j][1]*tmp[j][1])/NN/NN/s/s2;
					 tmp[j][1] = (
							 -out[j][0]*tmp[j][1]
							 +out[j][1]*tmpbuf)/NN/NN/s/s2;
				 }
				 // transform back into corr
				 fftw_execute(cr);
				 float maxcor=-1;
				 int poscor=0;
				 // use cross correlation for peak
				 for (int j=0; j < NN ; j++)
				 {
					 if (corr[j][0]>maxcor && (j < dvalue || NN-dvalue < j))
					 {
						 maxcor =corr[j][0];
						 poscor=j;
						 poscor=(j+NN/2)%NN-NN/2;

					 }
				 }
				 // we know maxcor now find the max in the indata
				 float maxin=-1.;
				 int maxpos =1;

				 for (int j=NN/2-dvalue; j < NN/2+dvalue ; j++)
				 {
					 if (in[j][0] > maxin)
					 {
						 maxpos = j; 
						 maxin = in[j][0];
					 }

				 }
				 for (int j=0; j < NN ; j++)
				 {
	//				if (Npeak%pvalue==17 && jvalue) 
//if (jvalue && abs(j-maxpos+200)<dvalue ) 
if (jvalue ) 
fprintf(corfile, "%8d %12.6f %12.6f %12.6f %d %d %d %d %12.6f %d\n", j-maxpos,in[j][0]/((maxin>0)?maxin:1),in2[j][0],corr[j][0],Npeak,shift,poscor,globalshift, maxin, maxpos);
// cat indata | plot 'u (int($5)%2==0?$1:NaN):2:5pal , "" u (int($5)%2==1?$1:NaN):(-$2):5 pal ; set xrange [-500:500]'

				 }
					 if (Npeak-lvalue==1) startshift = globalshift+maxpos;
					 if (Npeak>lvalue+1 && maxcor>cvalue && abs(maxcor)<dvalue) 
{
 fprintf(Npeak%2==tvalue?tickfile:tockfile,"%8d %5d %12.6f %d %d %d %d %d\n",Npeak,globalshift+(ovalue?maxpos:poscor)-startshift,maxcor,shift,poscor,maxpos,startshift,globalshift+poscor-startshift);
}
//cat tick | plot ' u 1:2  w lp pt 5 ps 2, "tock" u 1:2  w lp pt 5 ps 2 '


				 if (kvalue==0 && ((Npeak <10 || maxcor > 0.70) && Npeak%2==1))
				 {
					 shift = NN+poscor;//+globalshift/2;
				 }
				 else
				 {
					 shift=NN;
				 }
				 globalshift-=(NN-shift);
			 }
			 Npeak++;
		 }
		 fftw_destroy_plan(p);
		 fftw_destroy_plan(pr);
		 fftw_cleanup();
	 }
	 fclose(tickfile);
	 fclose(tockfile);
	 fclose(corfile);
	 fclose(rawfile);

     char command[1024] ;
     sprintf(command," echo 'uns colorbox; f=%d;h=%d/3600; set cbtics 1; set term png size 1920,1080 font \"DejaVuSansCondensed,12 truecolor \" ; set ytics nomirror; set out \"/dev/null\"; plot \"tick\" u ($1/h):($2/f); set y2tics ; set y2range [GPVAL_Y_MIN*f:GPVAL_Y_MAX*f]; set out \"/home/pi/lussen/www/tico.png\"; set fit quiet; set samples 1000; set ylabel sprintf(\"modulo 1/%%d s (s)\",h);\
d=15;d1=15; f(x)=(a-b*x+c*cos(3.1415926*(x-x0)*h/d)); g(x)=(a1-b1*x+c1*cos(3.1415926*(x-x1)*h/d1)); fit []f(x) \"tick\" u ($1/h):($2/f) via a,b,c,d,x0; print d,c*1000,\"ms \",b*86400,\"s/d\"; fit []g(x) \"tock\" u ($1/h):($2/f) via a1,b1,c1,d1,x1; print d1,c1*1000,\"ms \",b1*86400,\"s/d\"; set xrange [:]; set key below;set samples 1000; set xlabel \"time (s)\";\
set title sprintf(\"beaterror: %%.2fms\",1000*(a-a1));\
 plot \"tick\" u ($1/h):($2/f)   w p pt 13 ps 2  t sprintf(\"%%.1f s/d\",b*86400) , \"tock\" u ($1/h):($2/f)  w p pt 5 ps 2  t sprintf(\"%%.1f s/d\",b1*86400) , f(x) lc 1 t sprintf(\"ampl: %%.2fms\",c*1000), g(x) lc 2 t sprintf(\"ampl: %%.2fms\",c1*1000) ; print \"beaterror: \",1000*(a-a1), \"ms\" ; ' | gnuplot -persist ",fvalue,hvalue);


     if (vvalue) fprintf(stderr,"%s\n",command);
     return system(command);

}

/* functions */

struct HEADER readheader()
{
	int read = 0;
	unsigned char buffer4[4];
	unsigned char buffer2[2];

	struct HEADER header;
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
	
	 return header;
}


int openfiles(FILE **tickfile,  FILE **tockfile, 
              FILE **corfile,   FILE **rawfile, 
              FILE *pulsefile, FILE *pulsefiletock,
			  double *ps,      double *pst,
			  int jvalue,      int wvalue,
              int NN)
{
	int result = 0;
	float row;
	float val;
	*tickfile = fopen("tick", "w");
	if (*tickfile == 0)
	{
		fprintf (stderr,"cannot open tick\n");
		result++;
	}
	*tockfile = fopen("tock", "w");
	if (*tockfile == 0)
	{
		fprintf (stderr,"cannot open tock\n");
		result++;
	}
	if (jvalue)
	{
		*corfile = fopen("/media/verbext/indata", "w");
		if (*corfile == 0)
		{
			fprintf (stderr,"cannot open /media/verbext/indata\n");
			result++;
		}
	}
	if (wvalue)
	{
		*rawfile = fopen("raw", "w");
		if (*rawfile == 0)
		{
			fprintf (stderr,"cannot open raw\n");
			result++;
		}
	}
	pulsefile = fopen("pulseshape", "r");
	if (pulsefile == 0)
	{
		fprintf (stderr,"no pulsefile called pulseshape\n");
		result++;
	} 
	else
	{
		for (int p=0; p < NN ; p++) 
		{
			ps[p]=0.0;
			if (fscanf(pulsefile,"%g %g", &row,&val) != 2)
			{
				fprintf (stderr,"pulsefile should be 1 .0223 with %d rows\n",NN);
				result++;
			}
			ps[p] = val;
		}
		fclose(pulsefile);
	}
	pulsefiletock = fopen("pulseshapetock", "r");
	if (pulsefiletock == 0)
	{
		fprintf (stderr,"no pulsefile called pulseshapetock\n");
		result++;
	}
	else
	{
		for (int p=0; p < NN ; p++) 
		{
			pst[p]=0.0;
			if (fscanf(pulsefiletock,"%g %g", &row,&val) != 2)
			{
				fprintf (stderr,"pulsefiletock should be 1 .0223 with %d rows\n",NN);
				result++;
			}
			pst[p] = val;
		}
		fclose(pulsefiletock);
	}
	return result;
}
