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
              FILE *pulsefile, 
              FILE **tickavg, 
			  double *ps,      double *pst,
			  int jvalue, int wvalue, int svalue, int NN);

int main(int argc, char **argv) 
{
    FILE *tickavg=stdout;
    FILE *tickfile=stdout;
    FILE *tockfile=stdout;
    FILE *corfile=stdout;
    FILE *rawfile=stdout;
    FILE *pulsefile = 0 ;
    int nvalue = 48000;
    int dvalue = 4000;
    int c;
    int evalue = 0;
    float cvalue = 0;
    int ovalue = 0;
    int xvalue = -1;
    int pvalue = 15*2;
    int fvalue = 48000;
    int hvalue = 21600;
    int lvalue = 0;
    int rvalue = 0;
    int svalue = 0;
    int bvalue = 0;
    int uvalue = 0;
    int qvalue = 4000;
    int tvalue = 0;
    int kvalue = 0;
    int wvalue = 0;
    int vvalue = 0;
    int jvalue = 0;
    int NN     = 8000;
	int Ntick  = 0;
	int Ntock  = 0;
    opterr=0;
    double *ps;
    double *pst;
    double *avgtick;
    double *avgtock;
    int read = 0;
    opterr=0;

    while ((c = getopt (argc, argv, "n:d:l:r:q:twvh:f:e:op:jkc:sx:b:")) != -1)
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
            case 's':
                //split peaks
                svalue = 1;
                break;
            case 'j':
                //indata
                jvalue = 1;
                break;
            case 'b':
                // lower limit for filter
                bvalue = atoi(optarg);
                break;
            case 'x':
                // print one tiock
                xvalue = atoi(optarg);
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
                    fprintf (stderr, "h bph default 21600\nf wav frequency default 48000\nn maximum points\nd max distance for window shift\nl left trim (s)\nq move points up (default 2000)\nr right trim (s)\nj flatten the curve \nw raw input\nv show gnuplot command\ne gausfiliter stdev\np teethfor hisdev\nt toggle tick/tock\n s split tick and tock correlation peaks\nx <n> , print one tick/tock, completely\nb <n> bandpass all samples over <n>Hz\no use loudest noise and not correlation");
                }else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort ();
        }

    NN=(int)(fvalue*3600/hvalue);
    bvalue = (bvalue>0)?fvalue/bvalue:0;
	if (dvalue >NN/2) 
	{
		fprintf(stderr,"dvalue can be max half of the number of samples per tick\n");
		exit(-1);
	}
    pst = malloc(NN*sizeof(double));
    ps = malloc(NN*sizeof(double));
    avgtick = malloc(NN*sizeof(double));
    avgtock = malloc(NN*sizeof(double));

    lvalue = lvalue*hvalue/3600;
    rvalue = rvalue*hvalue/3600;
	if (openfiles(&tickfile,  &tockfile, 
				  &corfile,   &rawfile, 
				   pulsefile, &tickavg,
				   ps, pst, jvalue, wvalue,svalue,NN))
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


         p = fftw_plan_dft_1d(NN,  in,   out,       FFTW_FORWARD, FFTW_ESTIMATE );
         q = fftw_plan_dft_1d(NN,  in2,  filterFFT, FFTW_FORWARD, FFTW_ESTIMATE );
         pr = fftw_plan_dft_1d(NN, conv, in       , FFTW_BACKWARD, FFTW_ESTIMATE);
         cf = fftw_plan_dft_1d(NN, in2,  tmp,       FFTW_FORWARD, FFTW_ESTIMATE );
         cr = fftw_plan_dft_1d(NN, tmp,  corr,      FFTW_BACKWARD, FFTW_ESTIMATE);

		 if (evalue > 0)
         {
             // make filter array
             for (int j=0; j < NN; j++) 
             {
                 in2[j][0] = .398942280401/evalue*(exp(-((float)(j*j))/(float)(evalue*evalue)/2) + exp(-((float)(NN-j)*(NN-j))/(float)(evalue*evalue)/2));
                 in2[j][1] = 0.0;
             }

             fftw_execute(q);
         }
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
					 mean[j] = (float)((msb[0] << 8) | lsb[0]);
				 } else {
					 printf("Error reading file. %d bytes\n", read);
					 exit(-1);
				 }
			 }

			 for (int j=0; j < NN; j++) 
			 {

				 in[j][0] = (j>0)?abs(mean[j]-mean[j-1]):0.0;
				 in[j][1] = 0.0;
                 if (Npeak==xvalue) fprintf(rawfile, "%d %f %f\n",j,in[j][0],mean[j]);
			 }

			 // filter in array
			 if (evalue > 0)
			 {
				 fftw_execute(p);

				 for (int j=0; j < NN ; j++)
				 {
                     if (j<bvalue)
                     {
                         out[j][0]=0.0;
                         out[j][1]=0.0;
                     }

					 conv[j][0] = (
							 +out[j][0]*filterFFT[j][0]
							 -out[j][1]*filterFFT[j][1])/NN;
					 conv[j][1] = (
							 out[j][0]*filterFFT[j][1]
							 +out[j][1]*filterFFT[j][0])/NN;
				 }

				 fftw_execute(pr);
			 }

             if (Npeak==xvalue) for (int j=0; j < NN; j++) { fprintf(rawfile, "%d %f %f %f\n",j,in[j][0], filterFFT[j][0],filterFFT[j][1]); }
                     
			 float tot=0;
			 float mom=0;
			 int j=0;
			 int n=0;
			 if (Npeak>= lvalue)
			 {

				 double ix = 0.0;
				 double ixx =0.0;
				 for (int j=0; j < NN ; j++)
				 {
					 in[j][0] = (float)(in[j][0]/NN);
					 in[j][1] = 0.0;
					 ix+=in[j][0];
              //       if (wvalue) fprintf(rawfile, "%d %f\n",j,in[j][0]);
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
             if (Npeak==xvalue) for (int j=0; j < NN; j++) { fprintf(rawfile, "%d %f\n",j,corr[j][0]); }
				 float maxcor=-1;
				 int poscor=0;
				 // use cross correlation for peak
				 for (int j=0; j < NN ; j++)
				 {
					 if (corr[j][0]>maxcor && ( Npeak < lvalue+10 || (j < dvalue || NN-dvalue < j)))
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
             if (Npeak==xvalue) for (int j=0; j < NN; j++) { fprintf(rawfile, "%d %f\n",j-poscor,in[j][0]); }
				 for (int j=0; j < NN ; j++)
				 {
					 if (Npeak%2!=tvalue&&maxcor > cvalue) {
						 if (j-poscor >= 0 && j-poscor <NN)
						 {
							 if (j-poscor==0)Ntick++;
							 avgtick[j-poscor] += in[j][0];
						 }
					 }
					 if (Npeak%2==tvalue&&maxcor > cvalue) {
						 if (j-poscor >= 0 && j-poscor <NN)
						 {
							 if (j-poscor==0)Ntock++;
							 avgtock[j-poscor] += in[j][0];
						 }
					 }
	//				 				if (Npeak%pvalue==17 && jvalue) 
					 //if (jvalue && abs(j-maxpos+200)<dvalue ) 
			 if (jvalue ) 
						 fprintf(corfile, "%8d %12.6f %12.6f %12.6f %d %d %d %d %12.6f %d\n", j-(ovalue?maxpos:poscor),in[j][0]/((maxin>0)?maxin:1),in2[j][0],corr[j][0],Npeak,shift,poscor,globalshift, maxin, maxpos);
					 // cat indata | plot 'u (int($5)%2==0?$1:NaN):2:5pal , "" u (int($5)%2==1?$1:NaN):(-$2):5 pal ; set xrange [-500:500]'

				 }
					 if (Npeak-lvalue==1) startshift = globalshift+maxpos;
					 if (Npeak>lvalue+1 && maxcor>cvalue && abs(maxcor)<dvalue) 
				 {
 fprintf(Npeak%2==tvalue?tickfile:tockfile,"%8d %5d %12.6f %d %d %d %d %d\n",Npeak,globalshift+(ovalue?maxpos:poscor)-startshift,maxcor,shift,poscor,maxpos,startshift,globalshift+poscor-startshift);
				 }
//cat tick | plot ' u 1:2  w lp pt 5 ps 2, "tock" u 1:2  w lp pt 5 ps 2 '


				 if (kvalue==0 && ((Npeak < lvalue + 10 || maxcor > 0.70) && Npeak%2==1))
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
		for (int j=0; j < NN ; j++)
		{
			fprintf(tickavg,"%d %lf %lf \n",j,
				avgtick[j]/Ntick,
				avgtock[j]/Ntock);

//  cat shape  | plot 'u ($1-4000)/48.:2 w l t "tick" , "" u ($1-4000)/48.:3 w l t "tock" ; set xrange [-15:5]; set format y ""; set ylabel "abs(pressure)"; set xlabel "time (ms)"; set xtics 1 '

		}
	 }

	 fclose(tickavg);
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
     if (vvalue) fprintf(stderr," #cat shape  | plot 'u ($1-4000)/48.:($2*f($1)) w l   t \"tick snk795\" , \"\" u ($1-4000)/48.:($3*f($1)) w l t \"tock snk795\" ; set xrange [-14:2]; set ylabel \"abs(pressure) -mean\"; set xlabel \"time (ms)\"; set xtics 1 ; set x2tics ( \"200\" -14.058654,\"210\" -13.389194,\"220 \" -12.780594,\"230\" -12.224916,\"240\" -11.715545,\"250\" -11.246923,\"260\" -10.814349,\"270\" -10.413818 ,\"280\" -10.041896,\"290\" -9.695623,\"300\" -9.372436,\"310\" -9.070099,\"320\" -8.786659,\"330\" -8.520396 ,\"340\" -8.269796,\"350\" -8.033516,\"360\" -7.810363) rotate ;f(x)=(1+19*((x-4000)/48.<-8))'\n"); 
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
              FILE *pulsefile,
			  FILE **tickavg,
			  double *ps,      double *pst,
			  int jvalue,      int wvalue, int svalue,
              int NN)
{
	int result = 0;
	float row;
	float val,val2;

	*tickavg = fopen("shape", "w");
	if (*tickavg == 0)
	{
		fprintf (stderr,"cannot open shape\n");
		result++;
	}
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
int defaulttock[8000] = { -58, -79, -93, -115, -125, -140, -154, -171, -182, -181, -185, -187, -189, -188, -188, -188, -187, -187, -187, -187, -187, -187, -187, -187, -187, -187, -188, -188, -189, -189, -190, -191, -191, -192, -193, -194, -195, -196, -197, -198, -199, -200, -201, -202, -203, -204, -205, -206, -208, -209, -210, -211, -212, -213, -214, -215, -216, -217, -218, -218, -219, -220, -220, -221, -221, -221, -222, -222, -222, -222, -222, -222, -222, -222, -222, -222, -222, -222, -221, -221, -221, -221, -220, -220, -220, -220, -220, -219, -219, -219, -219, -219, -218, -218, -218, -218, -218, -218, -217, -217, -217, -217, -217, -216, -216, -216, -215, -215, -215, -214, -214, -213, -213, -212, -212, -211, -210, -210, -209, -208, -208, -207, -207, -206, -205, -205, -205, -204, -204, -204, -203, -203, -203, -203, -203, -204, -204, -205, -205, -206, -206, -207, -208, -209, -210, -211, -212, -213, -214, -215, -215, -216, -217, -218, -219, -219, -220, -220, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -220, -220, -220, -220, -219, -219, -219, -219, -219, -219, -219, -219, -219, -219, -219, -219, -219, -220, -220, -220, -221, -221, -222, -222, -223, -223, -224, -224, -225, -225, -225, -226, -226, -226, -226, -226, -227, -227, -227, -226, -226, -226, -226, -226, -226, -225, -225, -225, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -225, -225, -226, -226, -227, -228, -228, -229, -230, -231, -232, -233, -234, -235, -236, -237, -238, -239, -240, -241, -242, -243, -244, -245, -246, -247, -248, -249, -249, -250, -251, -251, -252, -252, -253, -253, -253, -253, -254, -254, -254, -254, -254, -253, -253, -253, -253, -252, -252, -252, -251, -251, -250, -250, -249, -249, -248, -247, -247, -246, -246, -245, -245, -244, -244, -243, -243, -242, -242, -241, -240, -240, -239, -239, -238, -238, -237, -236, -236, -235, -235, -234, -233, -232, -232, -231, -230, -230, -229, -228, -227, -227, -226, -225, -224, -224, -223, -222, -222, -221, -220, -220, -219, -219, -218, -218, -218, -218, -217, -217, -217, -217, -217, -218, -218, -218, -219, -219, -220, -221, -221, -222, -223, -224, -226, -227, -228, -229, -231, -232, -234, -235, -237, -238, -240, -242, -243, -245, -247, -248, -250, -251, -253, -254, -255, -257, -258, -259, -260, -261, -262, -263, -263, -264, -265, -265, -265, -265, -266, -266, -265, -265, -265, -265, -264, -264, -263, -262, -262, -261, -260, -259, -258, -257, -256, -255, -254, -253, -252, -251, -250, -249, -248, -247, -246, -245, -244, -243, -242, -242, -241, -241, -240, -240, -240, -240, -240, -240, -240, -240, -241, -241, -242, -242, -243, -244, -245, -246, -247, -248, -249, -250, -251, -253, -254, -255, -256, -258, -259, -260, -261, -262, -263, -264, -265, -266, -266, -267, -268, -268, -268, -269, -269, -269, -269, -269, -268, -268, -268, -267, -267, -266, -265, -264, -264, -263, -262, -261, -260, -259, -258, -257, -256, -255, -255, -254, -253, -252, -252, -251, -251, -251, -250, -250, -250, -250, -250, -250, -250, -251, -251, -251, -252, -252, -253, -253, -254, -254, -255, -256, -256, -257, -257, -258, -258, -258, -259, -259, -259, -259, -259, -259, -259, -259, -259, -258, -258, -258, -257, -256, -256, -255, -254, -253, -253, -252, -251, -250, -249, -248, -247, -246, -245, -245, -244, -243, -242, -241, -240, -239, -238, -237, -236, -236, -235, -234, -233, -232, -231, -230, -229, -228, -227, -226, -225, -225, -224, -223, -222, -222, -221, -221, -220, -220, -219, -219, -219, -219, -219, -220, -220, -221, -221, -222, -223, -224, -225, -226, -227, -228, -229, -231, -232, -234, -235, -237, -238, -240, -242, -243, -245, -247, -248, -250, -251, -253, -254, -256, -257, -259, -260, -261, -262, -263, -264, -265, -266, -267, -267, -268, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -268, -268, -268, -267, -266, -266, -265, -265, -264, -263, -263, -262, -261, -261, -260, -260, -259, -259, -259, -259, -258, -258, -258, -258, -258, -258, -259, -259, -259, -259, -260, -260, -260, -261, -261, -261, -262, -262, -262, -262, -262, -262, -262, -262, -262, -261, -261, -261, -260, -260, -259, -258, -258, -257, -256, -255, -254, -253, -252, -251, -250, -249, -248, -247, -246, -245, -245, -244, -243, -242, -242, -241, -241, -240, -240, -240, -239, -239, -239, -239, -239, -240, -240, -240, -241, -242, -242, -243, -244, -245, -246, -247, -248, -249, -250, -252, -253, -254, -255, -257, -258, -259, -261, -262, -263, -264, -266, -267, -268, -269, -270, -271, -272, -272, -273, -273, -274, -274, -275, -275, -275, -275, -275, -275, -275, -275, -274, -274, -274, -273, -273, -272, -272, -271, -270, -270, -269, -268, -267, -267, -266, -265, -265, -264, -264, -263, -263, -263, -262, -262, -262, -262, -262, -262, -262, -262, -262, -262, -263, -263, -263, -264, -264, -265, -265, -266, -266, -266, -267, -267, -268, -268, -268, -268, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -269, -270, -270, -271, -271, -272, -272, -273, -274, -275, -275, -276, -277, -278, -279, -280, -281, -282, -283, -284, -285, -285, -286, -287, -287, -288, -288, -288, -289, -289, -289, -289, -288, -288, -288, -287, -286, -286, -285, -284, -283, -282, -281, -280, -279, -277, -276, -275, -273, -272, -271, -269, -268, -266, -265, -264, -263, -261, -260, -259, -258, -257, -256, -255, -255, -254, -254, -253, -253, -253, -252, -252, -252, -252, -252, -252, -253, -253, -253, -253, -254, -254, -254, -255, -255, -256, -256, -256, -257, -257, -257, -258, -258, -258, -259, -259, -259, -260, -260, -260, -261, -261, -262, -262, -262, -263, -263, -264, -265, -265, -266, -267, -267, -268, -269, -269, -270, -271, -272, -272, -273, -274, -274, -275, -276, -276, -277, -277, -277, -277, -277, -277, -277, -277, -277, -276, -276, -275, -275, -274, -273, -272, -271, -270, -269, -267, -266, -265, -264, -262, -261, -260, -258, -257, -256, -254, -253, -252, -251, -250, -250, -249, -248, -248, -247, -247, -247, -247, -246, -246, -247, -247, -247, -247, -248, -248, -249, -249, -250, -250, -251, -252, -252, -253, -253, -254, -255, -255, -256, -256, -257, -257, -258, -258, -258, -259, -259, -259, -260, -260, -260, -261, -261, -261, -262, -262, -263, -263, -263, -264, -264, -265, -266, -266, -267, -268, -268, -269, -270, -271, -272, -273, -275, -276, -277, -278, -280, -281, -283, -284, -286, -287, -289, -291, -292, -294, -296, -297, -299, -301, -302, -304, -306, -307, -309, -310, -312, -313, -314, -316, -317, -318, -319, -320, -321, -322, -322, -323, -323, -324, -324, -324, -325, -325, -325, -325, -325, -325, -324, -324, -324, -324, -324, -323, -323, -323, -322, -322, -321, -321, -321, -320, -320, -319, -319, -318, -318, -317, -316, -316, -315, -314, -313, -312, -311, -310, -309, -308, -307, -306, -304, -303, -302, -300, -299, -298, -296, -295, -294, -292, -291, -290, -289, -288, -287, -286, -286, -285, -285, -284, -284, -284, -284, -284, -284, -284, -284, -284, -285, -285, -285, -286, -286, -286, -287, -287, -287, -287, -287, -287, -287, -287, -286, -286, -286, -285, -285, -284, -284, -283, -282, -282, -281, -281, -280, -279, -279, -279, -278, -278, -278, -277, -277, -277, -277, -277, -278, -278, -278, -278, -279, -279, -279, -280, -280, -281, -281, -281, -282, -282, -282, -282, -283, -283, -283, -283, -283, -283, -282, -282, -282, -282, -281, -281, -280, -280, -279, -279, -278, -278, -277, -277, -276, -276, -275, -275, -275, -274, -274, -273, -273, -273, -273, -273, -273, -272, -272, -272, -273, -273, -273, -273, -273, -274, -274, -274, -275, -275, -276, -276, -277, -277, -278, -278, -279, -280, -280, -281, -282, -282, -283, -284, -284, -285, -286, -286, -287, -288, -288, -289, -290, -290, -291, -292, -292, -293, -294, -294, -295, -295, -296, -297, -297, -298, -299, -299, -300, -301, -301, -302, -302, -303, -304, -304, -305, -305, -306, -306, -306, -307, -307, -307, -307, -307, -307, -307, -307, -307, -307, -306, -306, -305, -305, -304, -303, -303, -302, -301, -300, -299, -298, -297, -296, -295, -294, -293, -292, -291, -290, -289, -288, -288, -287, -286, -285, -285, -284, -284, -283, -283, -283, -282, -282, -282, -282, -282, -282, -282, -282, -282, -282, -283, -283, -283, -283, -284, -284, -285, -285, -285, -286, -286, -287, -287, -288, -288, -289, -289, -290, -290, -291, -291, -292, -292, -293, -293, -293, -294, -294, -295, -295, -296, -296, -296, -297, -297, -297, -297, -298, -298, -298, -298, -298, -298, -298, -298, -298, -297, -297, -297, -297, -296, -296, -295, -295, -295, -294, -294, -294, -293, -293, -293, -293, -292, -292, -292, -292, -292, -292, -292, -292, -293, -293, -293, -293, -294, -294, -294, -295, -295, -296, -296, -297, -297, -298, -298, -299, -299, -300, -301, -301, -302, -302, -303, -303, -303, -304, -304, -305, -305, -305, -306, -306, -306, -306, -306, -306, -306, -306, -306, -306, -306, -306, -306, -306, -305, -305, -305, -304, -304, -304, -303, -303, -303, -302, -302, -302, -301, -301, -301, -300, -300, -300, -299, -299, -299, -299, -299, -299, -299, -298, -298, -298, -299, -299, -299, -299, -299, -299, -299, -299, -299, -300, -300, -300, -300, -300, -300, -300, -300, -300, -300, -300, -300, -299, -299, -299, -299, -299, -299, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -298, -297, -297, -297, -297, -297, -297, -296, -296, -295, -295, -294, -294, -293, -292, -291, -290, -290, -289, -288, -287, -286, -285, -284, -283, -282, -281, -280, -279, -278, -278, -277, -276, -276, -275, -275, -275, -275, -274, -274, -274, -274, -275, -275, -275, -276, -276, -276, -277, -278, -278, -279, -280, -280, -281, -282, -283, -283, -284, -285, -286, -286, -287, -288, -288, -289, -290, -290, -291, -291, -292, -292, -292, -293, -293, -294, -294, -294, -295, -295, -295, -296, -296, -297, -297, -298, -298, -299, -300, -300, -301, -302, -303, -304, -305, -306, -307, -308, -309, -310, -311, -312, -313, -314, -315, -315, -316, -317, -318, -319, -319, -320, -320, -321, -321, -322, -322, -322, -322, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -323, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -324, -323, -323, -323, -323, -323, -322, -322, -322, -322, -321, -321, -321, -320, -320, -320, -319, -319, -319, -318, -318, -317, -317, -317, -317, -316, -316, -316, -315, -315, -315, -315, -315, -315, -314, -314, -314, -314, -314, -314, -313, -313, -313, -313, -313, -312, -312, -312, -311, -311, -310, -310, -309, -309, -308, -308, -307, -307, -306, -306, -305, -305, -304, -304, -304, -304, -303, -303, -303, -304, -304, -304, -305, -305, -306, -306, -307, -308, -309, -310, -311, -312, -313, -314, -315, -316, -317, -318, -319, -320, -321, -322, -323, -324, -325, -326, -326, -327, -328, -328, -329, -329, -330, -330, -331, -331, -332, -332, -332, -333, -333, -333, -334, -334, -335, -335, -336, -336, -337, -337, -338, -338, -339, -340, -340, -341, -342, -342, -343, -343, -344, -344, -345, -345, -346, -346, -346, -346, -346, -346, -346, -346, -345, -344, -344, -343, -342, -341, -340, -338, -337, -335, -334, -332, -330, -328, -327, -325, -323, -320, -318, -316, -314, -312, -310, -308, -306, -304, -302, -300, -298, -297, -295, -293, -292, -290, -289, -288, -287, -286, -285, -284, -283, -282, -281, -281, -280, -280, -280, -279, -279, -279, -279, -279, -279, -279, -279, -279, -279, -279, -279, -279, -280, -280, -280, -281, -281, -282, -282, -283, -283, -284, -284, -285, -286, -286, -287, -288, -288, -289, -290, -291, -291, -292, -293, -294, -294, -295, -295, -296, -297, -297, -297, -298, -298, -298, -298, -299, -299, -299, -298, -298, -298, -297, -297, -297, -296, -295, -295, -294, -293, -292, -291, -290, -289, -288, -287, -286, -285, -284, -283, -282, -282, -281, -280, -279, -279, -278, -277, -277, -276, -276, -276, -276, -276, -276, -276, -276, -276, -276, -277, -277, -278, -278, -279, -280, -280, -281, -282, -283, -284, -285, -286, -287, -288, -289, -290, -291, -292, -293, -294, -295, -296, -297, -298, -299, -300, -301, -301, -302, -302, -303, -303, -304, -304, -304, -304, -304, -304, -304, -304, -304, -304, -303, -303, -303, -303, -302, -302, -302, -302, -302, -302, -302, -302, -302, -303, -303, -303, -304, -304, -305, -306, -307, -307, -308, -309, -310, -311, -313, -314, -315, -316, -317, -318, -319, -321, -322, -323, -324, -325, -326, -327, -327, -328, -329, -329, -330, -330, -330, -331, -331, -331, -331, -330, -330, -330, -329, -329, -328, -328, -327, -326, -326, -325, -324, -323, -323, -322, -321, -320, -320, -319, -319, -318, -318, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -317, -318, -318, -318, -318, -318, -318, -318, -317, -317, -317, -317, -317, -317, -317, -316, -316, -316, -316, -316, -315, -315, -315, -315, -314, -314, -314, -313, -313, -313, -312, -312, -312, -311, -311, -311, -310, -310, -310, -309, -309, -308, -308, -307, -307, -306, -306, -305, -305, -305, -304, -304, -303, -303, -303, -302, -302, -302, -301, -301, -301, -301, -300, -300, -300, -300, -300, -300, -300, -300, -300, -301, -301, -301, -301, -302, -302, -302, -303, -303, -304, -305, -305, -306, -307, -308, -308, -309, -310, -311, -312, -313, -314, -315, -316, -317, -318, -319, -320, -321, -321, -322, -323, -324, -325, -326, -326, -327, -328, -329, -329, -330, -330, -331, -331, -331, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -332, -331, -331, -331, -331, -331, -331, -331, -330, -330, -330, -330, -330, -331, -331, -331, -331, -332, -332, -332, -333, -334, -334, -335, -335, -336, -337, -338, -339, -339, -340, -341, -342, -343, -344, -345, -346, -347, -348, -349, -350, -351, -352, -353, -354, -355, -356, -357, -358, -359, -360, -361, -362, -362, -363, -364, -364, -365, -365, -366, -366, -367, -367, -367, -367, -367, -367, -367, -367, -366, -366, -365, -365, -364, -363, -362, -361, -360, -359, -358, -356, -355, -353, -351, -350, -348, -346, -344, -342, -340, -338, -336, -334, -332, -329, -327, -325, -323, -321, -319, -317, -315, -313, -312, -310, -309, -307, -306, -305, -304, -303, -303, -302, -302, -301, -301, -301, -302, -302, -303, -303, -304, -305, -306, -307, -308, -310, -311, -312, -314, -316, -317, -319, -320, -322, -324, -326, -327, -329, -331, -332, -334, -335, -337, -338, -339, -341, -342, -343, -344, -345, -346, -347, -348, -348, -349, -349, -350, -350, -351, -351, -351, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -351, -351, -351, -351, -350, -350, -350, -349, -349, -348, -348, -347, -346, -346, -345, -344, -344, -343, -342, -341, -340, -339, -338, -338, -337, -336, -335, -334, -333, -332, -331, -330, -329, -328, -327, -326, -325, -324, -322, -321, -320, -320, -319, -318, -317, -316, -315, -314, -313, -313, -312, -311, -311, -310, -310, -309, -309, -309, -309, -309, -309, -309, -310, -310, -310, -311, -312, -312, -313, -314, -315, -316, -317, -318, -320, -321, -322, -323, -324, -326, -327, -328, -329, -330, -331, -332, -333, -333, -334, -335, -335, -335, -336, -336, -336, -336, -336, -335, -335, -335, -334, -334, -333, -333, -332, -331, -331, -330, -329, -329, -328, -328, -327, -327, -327, -326, -326, -326, -326, -326, -326, -326, -327, -327, -327, -328, -328, -329, -329, -330, -331, -331, -332, -333, -333, -334, -334, -335, -335, -336, -336, -336, -337, -337, -337, -337, -337, -336, -336, -336, -336, -335, -335, -334, -334, -333, -333, -332, -331, -331, -330, -330, -329, -329, -329, -328, -328, -328, -328, -328, -328, -328, -328, -329, -329, -329, -330, -331, -331, -332, -333, -334, -335, -336, -337, -338, -339, -339, -340, -341, -342, -342, -343, -343, -343, -343, -344, -343, -343, -343, -342, -342, -341, -340, -339, -338, -337, -336, -335, -334, -332, -331, -330, -329, -327, -326, -325, -324, -322, -321, -320, -319, -318, -318, -317, -316, -316, -315, -315, -314, -314, -314, -314, -314, -314, -314, -314, -314, -315, -315, -316, -316, -317, -317, -318, -319, -320, -321, -322, -323, -324, -325, -326, -327, -328, -330, -331, -332, -333, -335, -336, -337, -339, -340, -341, -342, -343, -344, -345, -346, -347, -348, -349, -349, -350, -350, -351, -351, -351, -351, -352, -352, -351, -351, -351, -351, -351, -351, -351, -350, -350, -350, -350, -350, -350, -350, -350, -350, -350, -350, -351, -351, -352, -352, -353, -353, -354, -355, -355, -356, -357, -358, -358, -359, -360, -361, -361, -362, -363, -364, -364, -365, -365, -366, -366, -367, -367, -368, -368, -368, -369, -369, -369, -369, -369, -370, -370, -370, -370, -370, -370, -370, -369, -369, -369, -369, -369, -368, -368, -367, -367, -366, -366, -365, -364, -364, -363, -362, -361, -360, -359, -358, -357, -355, -354, -353, -352, -350, -349, -348, -347, -345, -344, -343, -342, -340, -339, -338, -337, -336, -335, -335, -334, -333, -333, -332, -332, -332, -332, -332, -332, -332, -332, -332, -333, -333, -334, -334, -335, -336, -336, -337, -338, -339, -339, -340, -341, -342, -343, -344, -344, -345, -346, -347, -347, -348, -349, -349, -350, -350, -351, -351, -351, -352, -352, -352, -352, -352, -353, -353, -352, -352, -352, -352, -352, -351, -351, -350, -350, -349, -348, -348, -347, -346, -345, -344, -344, -343, -342, -341, -340, -339, -338, -337, -337, -336, -335, -334, -333, -333, -332, -331, -331, -330, -330, -329, -329, -328, -328, -327, -327, -326, -326, -325, -324, -324, -323, -323, -322, -322, -321, -320, -320, -319, -319, -318, -318, -317, -317, -316, -316, -316, -316, -316, -316, -316, -316, -316, -316, -317, -317, -318, -319, -319, -320, -321, -322, -323, -325, -326, -327, -328, -330, -331, -332, -334, -335, -336, -338, -339, -340, -342, -343, -344, -345, -346, -347, -348, -349, -349, -350, -351, -351, -351, -352, -352, -352, -352, -352, -352, -352, -352, -352, -351, -351, -351, -350, -350, -349, -349, -349, -348, -348, -348, -347, -347, -347, -347, -346, -346, -346, -346, -346, -346, -346, -347, -347, -347, -347, -348, -348, -348, -349, -349, -350, -350, -350, -351, -351, -351, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -352, -351, -351, -351, -351, -350, -350, -350, -350, -349, -349, -349, -349, -348, -348, -348, -348, -348, -348, -348, -347, -347, -347, -347, -347, -347, -347, -347, -347, -347, -347, -346, -346, -346, -346, -346, -346, -345, -345, -345, -345, -344, -344, -344, -344, -343, -343, -343, -343, -342, -342, -342, -342, -342, -341, -341, -341, -341, -341, -341, -341, -340, -340, -340, -340, -340, -340, -340, -340, -340, -340, -339, -339, -339, -339, -339, -338, -338, -338, -337, -337, -337, -336, -336, -335, -334, -334, -333, -333, -332, -332, -331, -330, -330, -329, -329, -329, -328, -328, -328, -327, -327, -327, -327, -327, -327, -327, -328, -328, -328, -328, -329, -329, -330, -330, -331, -331, -332, -333, -333, -334, -335, -335, -336, -337, -337, -338, -339, -339, -340, -341, -341, -342, -343, -344, -344, -345, -346, -346, -347, -348, -349, -349, -350, -351, -352, -352, -353, -354, -354, -355, -356, -357, -357, -358, -358, -359, -360, -360, -361, -361, -362, -362, -363, -363, -363, -364, -364, -364, -364, -364, -364, -364, -364, -364, -364, -364, -364, -363, -363, -362, -362, -361, -360, -360, -359, -358, -357, -356, -355, -354, -353, -352, -350, -349, -348, -346, -345, -343, -342, -341, -339, -337, -336, -334, -333, -331, -330, -328, -327, -325, -324, -323, -321, -320, -319, -318, -317, -316, -315, -314, -313, -313, -312, -312, -311, -311, -311, -311, -311, -311, -312, -312, -312, -313, -314, -314, -315, -316, -317, -318, -319, -320, -322, -323, -324, -326, -327, -328, -330, -331, -333, -334, -336, -337, -339, -340, -342, -343, -345, -347, -348, -350, -351, -353, -354, -356, -357, -359, -360, -361, -363, -364, -365, -366, -367, -368, -369, -370, -371, -371, -372, -372, -373, -373, -373, -374, -374, -374, -374, -373, -373, -373, -373, -372, -372, -371, -371, -370, -369, -369, -368, -367, -367, -366, -365, -364, -363, -362, -361, -361, -360, -359, -358, -357, -356, -356, -355, -354, -354, -353, -353, -352, -352, -351, -351, -351, -351, -351, -351, -351, -351, -351, -352, -352, -353, -353, -354, -354, -355, -356, -356, -357, -358, -359, -360, -361, -361, -362, -363, -364, -365, -365, -366, -367, -367, -368, -368, -369, -369, -370, -370, -370, -371, -371, -371, -371, -371, -372, -372, -372, -372, -372, -373, -373, -373, -374, -374, -375, -375, -376, -377, -377, -378, -379, -380, -381, -382, -383, -384, -386, -387, -388, -389, -391, -392, -394, -395, -397, -398, -399, -401, -402, -404, -405, -406, -407, -408, -409, -410, -411, -412, -413, -413, -414, -414, -414, -414, -414, -414, -413, -412, -412, -411, -409, -408, -406, -405, -403, -400, -398, -395, -392, -389, -386, -382, -378, -374, -370, -365, -360, -354, -349, -343, -336, -329, -322, -315, -307, -299, -290, -281, -271, -262, -251, -241, -230, -219, -207, -195, -183, -171, -158, -146, -133, -120, -107, -94, -81, -69, -56, -44, -32, -20, -8, 3, 14, 24, 34, 43, 52, 60, 68, 75, 82, 88, 93, 98, 102, 106, 109, 112, 114, 116, 117, 118, 119, 120, 120, 120, 120, 120, 120, 120, 120, 121, 121, 121, 122, 123, 124, 126, 127, 129, 132, 134, 137, 140, 143, 147, 151, 155, 159, 163, 167, 172, 176, 181, 185, 189, 193, 197, 201, 205, 208, 211, 214, 216, 218, 220, 221, 222, 223, 223, 223, 222, 221, 219, 218, 215, 213, 210, 206, 202, 198, 194, 189, 184, 179, 173, 167, 161, 154, 147, 140, 133, 126, 118, 110, 102, 94, 86, 78, 70, 61, 53, 44, 36, 27, 18, 10, 1, -7, -15, -24, -32, -40, -48, -56, -64, -71, -79, -86, -93, -100, -107, -114, -121, -127, -133, -139, -145, -151, -156, -162, -167, -172, -176, -181, -185, -189, -193, -196, -199, -202, -205, -207, -209, -210, -212, -212, -213, -213, -212, -211, -210, -208, -205, -202, -198, -193, -188, -182, -175, -167, -158, -148, -137, -124, -111, -96, -79, -61, -41, -19, 5, 31, 59, 90, 124, 161, 200, 244, 290, 341, 396, 455, 519, 588, 662, 742, 828, 920, 1018, 1124, 1237, 1357, 1485, 1621, 1764, 1916, 2077, 2246, 2423, 2608, 2801, 3002, 3210, 3426, 3647, 3875, 4108, 4345, 4585, 4828, 5072, 5317, 5560, 5801, 6039, 6272, 6499, 6718, 6928, 7129, 7318, 7495, 7659, 7808, 7941, 8058, 8159, 8242, 8308, 8355, 8385, 8396, 8390, 8366, 8325, 8267, 8193, 8104, 8001, 7884, 7754, 7613, 7462, 7301, 7132, 6955, 6773, 6585, 6393, 6198, 6002, 5803, 5605, 5407, 5210, 5015, 4822, 4632, 4446, 4264, 4086, 3913, 3745, 3583, 3425, 3274, 3128, 2988, 2854, 2726, 2604, 2488, 2378, 2274, 2176, 2083, 1997, 1915, 1839, 1768, 1702, 1641, 1585, 1533, 1485, 1440, 1400, 1362, 1328, 1297, 1268, 1241, 1216, 1193, 1172, 1151, 1132, 1113, 1095, 1078, 1061, 1043, 1026, 1009, 991, 974, 956, 938, 919, 901, 882, 863, 845, 826, 807, 789, 771, 754, 737, 721, 706, 691, 677, 665, 653, 643, 633, 625, 618, 612, 607, 603, 600, 597, 596, 595, 595, 596, 597, 598, 599, 600, 602, 603, 604, 605, 605, 605, 605, 603, 602, 599, 596, 593, 588, 583, 578, 572, 565, 558, 551, 543, 535, 527, 519, 510, 502, 493, 485, 476, 468, 460, 452, 444, 436, 429, 422, 415, 408, 402, 395, 389, 383, 377, 371, 365, 359, 353, 347, 340, 334, 327, 320, 313, 306, 298, 291, 283, 275, 266, 258, 249, 240, 231, 221, 212, 203, 193, 184, 175, 166, 157, 148, 140, 131, 123, 116, 109, 102, 95, 89, 84, 79, 74, 70, 66, 63, 60, 58, 57, 55, 55, 55, 55, 56, 57, 59, 61, 64, 68, 72, 76, 82, 88, 94, 102, 110, 120, 130, 142, 155, 170, 186, 204, 224, 246, 271, 298, 328, 362, 398, 439, 483, 531, 583, 640, 702, 769, 840, 918, 1001, 1089, 1183, 1283, 1389, 1500, 1617, 1740, 1868, 2001, 2139, 2282, 2428, 2579, 2734, 2891, 3050, 3211, 3374, 3537, 3700, 3862, 4023, 4182, 4338, 4490, 4639, 4782, 4921, 5053, 5179, 5298, 5409, 5513, 5609, 5696, 5774, 5844, 5905, 5957, 6000, 6035, 6060, 6077, 6086, 6087, 6081, 6067, 6046, 6018, 5985, 5946, 5901, 5852, 5799, 5741, 5680, 5617, 5550, 5481, 5411, 5338, 5265, 5190, 5115, 5039, 4963, 4887, 4811, 4735, 4660, 4585, 4511, 4437, 4364, 4292, 4220, 4150, 4080, 4012, 3944, 3877, 3811, 3747, 3683, 3621, 3560, 3500, 3441, 3383, 3327, 3272, 3218, 3166, 3115, 3065, 3017, 2970, 2925, 2881, 2838, 2796, 2756, 2717, 2679, 2643, 2607, 2572, 2538, 2505, 2472, 2440, 2409, 2378, 2347, 2316, 2286, 2255, 2225, 2194, 2163, 2131, 2100, 2068, 2035, 2003, 1969, 1936, 1902, 1867, 1833, 1798, 1762, 1727, 1691, 1655, 1620, 1584, 1548, 1513, 1478, 1443, 1409, 1375, 1341, 1309, 1277, 1246, 1215, 1185, 1157, 1129, 1102, 1076, 1051, 1027, 1004, 983, 962, 942, 923, 906, 889, 873, 858, 845, 832, 820, 808, 798, 788, 779, 771, 763, 756, 750, 743, 738, 732, 727, 723, 718, 714, 709, 705, 701, 697, 692, 688, 683, 678, 674, 668, 663, 657, 651, 645, 639, 632, 625, 618, 611, 603, 595, 587, 579, 570, 561, 553, 544, 535, 526, 517, 508, 498, 489, 480, 471, 462, 453, 444, 435, 426, 417, 409, 400, 392, 384, 376, 368, 361, 353, 346, 339, 333, 326, 320, 314, 308, 303, 298, 293, 289, 284, 280, 276, 272, 269, 266, 263, 260, 257, 254, 252, 249, 247, 245, 242, 240, 238, 236, 234, 232, 229, 227, 225, 223, 221, 218, 216, 214, 211, 209, 207, 204, 202, 199, 197, 194, 192, 189, 187, 184, 182, 179, 177, 175, 172, 170, 168, 166, 163, 161, 159, 157, 156, 154, 152, 150, 149, 147, 146, 144, 143, 142, 140, 139, 138, 137, 136, 135, 135, 134, 133, 133, 132, 132, 131, 131, 130, 130, 130, 129, 129, 129, 128, 128, 128, 127, 127, 126, 126, 125, 125, 124, 123, 123, 122, 121, 120, 119, 118, 118, 117, 116, 115, 114, 113, 112, 111, 110, 109, 108, 108, 107, 106, 105, 105, 104, 103, 103, 102, 102, 101, 100, 100, 100, 99, 99, 98, 98, 98, 97, 97, 97, 97, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 100, 100, 100, 99, 99, 99, 99, 99, 98, 98, 98, 97, 97, 96, 96, 95, 95, 94, 94, 93, 93, 92, 92, 91, 91, 91, 90, 90, 90, 90, 90, 90, 91, 91, 91, 92, 93, 94, 95, 96, 97, 98, 99, 101, 102, 104, 105, 107, 109, 110, 112, 113, 115, 117, 118, 119, 121, 122, 123, 125, 126, 127, 128, 129, 130, 130, 131, 132, 132, 133, 133, 134, 134, 135, 135, 135, 135, 135, 136, 136, 136, 135, 135, 135, 135, 135, 134, 134, 133, 133, 132, 131, 131, 130, 129, 129, 128, 127, 126, 125, 124, 124, 123, 122, 121, 120, 119, 119, 118, 117, 116, 116, 115, 115, 114, 114, 114, 113, 113, 113, 113, 114, 114, 114, 115, 116, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 135, 136, 136, 137, 137, 137, 136, 136, 136, 135, 135, 134, 133, 133, 132, 131, 130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 121, 120, 119, 118, 118, 117, 117, 116, 116, 116, 115, 115, 115, 115, 115, 115, 114, 114, 114, 114, 114, 114, 115, 115, 115, 115, 115, 115, 115, 115, 116, 116, 116, 116, 117, 117, 118, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 129, 130, 131, 133, 134, 136, 137, 139, 140, 142, 143, 144, 146, 147, 148, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 184, 185, 186, 187, 188, 188, 189, 190, 190, 191, 192, 192, 193, 193, 193, 194, 194, 194, 195, 195, 195, 195, 195, 194, 194, 194, 193, 193, 192, 192, 191, 190, 189, 188, 187, 186, 185, 184, 183, 181, 180, 179, 178, 176, 175, 173, 172, 171, 169, 168, 167, 165, 164, 163, 161, 160, 159, 158, 157, 156, 155, 154, 153, 152, 152, 151, 150, 150, 149, 149, 149, 148, 148, 148, 148, 148, 147, 147, 147, 148, 148, 148, 148, 148, 148, 149, 149, 149, 150, 150, 150, 151, 151, 152, 153, 153, 154, 155, 155, 156, 157, 158, 159, 159, 160, 161, 162, 163, 164, 165, 166, 168, 169, 170, 171, 172, 173, 175, 176, 177, 178, 180, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 196, 197, 197, 198, 198, 199, 199, 199, 199, 200, 200, 200, 200, 200, 200, 200, 200, 200, 201, 201, 201, 201, 201, 202, 202, 202, 203, 203, 204, 204, 204, 205, 205, 206, 206, 207, 207, 207, 208, 208, 208, 208, 209, 209, 209, 209, 208, 208, 208, 208, 207, 207, 206, 206, 205, 204, 204, 203, 202, 201, 200, 199, 199, 198, 197, 196, 195, 194, 194, 193, 193, 192, 192, 191, 191, 191, 190, 190, 190, 190, 190, 190, 190, 190, 191, 191, 191, 191, 191, 190, 190, 190, 190, 189, 189, 188, 187, 186, 186, 184, 183, 182, 181, 179, 178, 176, 174, 173, 171, 169, 167, 165, 163, 161, 158, 156, 154, 152, 150, 148, 146, 144, 142, 140, 139, 137, 136, 134, 133, 132, 131, 130, 129, 129, 129, 128, 128, 128, 129, 129, 130, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 143, 144, 145, 147, 148, 150, 151, 153, 154, 156, 157, 159, 160, 162, 163, 164, 166, 167, 168, 170, 171, 172, 173, 174, 175, 176, 177, 178, 178, 179, 179, 180, 180, 180, 180, 180, 179, 179, 179, 178, 177, 177, 176, 175, 174, 173, 172, 171, 170, 168, 167, 166, 165, 164, 163, 162, 161, 161, 160, 160, 159, 159, 159, 158, 158, 158, 159, 159, 159, 160, 160, 161, 161, 162, 163, 163, 164, 164, 165, 166, 166, 167, 167, 167, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 169, 169, 169, 170, 170, 170, 171, 171, 172, 173, 173, 174, 174, 175, 175, 176, 176, 177, 177, 178, 178, 178, 179, 179, 179, 179, 179, 178, 178, 178, 177, 176, 176, 175, 174, 173, 172, 171, 170, 169, 168, 166, 165, 164, 163, 162, 161, 159, 158, 157, 156, 156, 155, 154, 153, 153, 152, 152, 152, 151, 151, 151, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 151, 151, 151, 151, 151, 152, 152, 152, 152, 153, 153, 153, 154, 154, 154, 154, 155, 155, 155, 155, 155, 155, 156, 156, 156, 156, 155, 155, 155, 155, 155, 155, 155, 154, 154, 154, 153, 153, 153, 152, 152, 151, 151, 150, 150, 149, 148, 147, 147, 146, 145, 144, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 132, 131, 130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 119, 118, 117, 117, 116, 116, 116, 115, 115, 115, 115, 115, 116, 116, 116, 117, 118, 119, 119, 120, 122, 123, 124, 125, 126, 128, 129, 131, 132, 133, 135, 136, 137, 138, 140, 141, 142, 142, 143, 144, 144, 145, 145, 145, 145, 145, 144, 144, 143, 142, 141, 140, 138, 137, 135, 134, 132, 130, 128, 126, 124, 122, 120, 118, 116, 114, 112, 110, 108, 106, 104, 102, 101, 99, 98, 96, 95, 94, 93, 92, 91, 90, 89, 88, 88, 87, 87, 86, 86, 85, 85, 85, 85, 84, 84, 84, 84, 84, 84, 85, 85, 85, 85, 86, 86, 86, 87, 87, 88, 88, 88, 89, 89, 89, 90, 90, 90, 90, 90, 90, 90, 89, 89, 88, 88, 87, 86, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 73, 72, 71, 70, 69, 69, 68, 68, 67, 67, 67, 67, 66, 66, 66, 67, 67, 67, 67, 68, 68, 69, 69, 70, 71, 71, 72, 73, 74, 75, 75, 76, 77, 78, 79, 79, 80, 81, 81, 82, 83, 83, 84, 84, 85, 85, 86, 86, 87, 87, 88, 89, 89, 90, 90, 91, 91, 92, 92, 93, 93, 94, 94, 95, 95, 96, 96, 97, 97, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99, 99, 98, 98, 97, 97, 96, 95, 94, 93, 92, 91, 90, 88, 87, 86, 84, 82, 81, 79, 78, 76, 74, 73, 71, 70, 68, 67, 66, 65, 64, 63, 62, 62, 61, 61, 60, 60, 60, 60, 60, 60, 61, 61, 61, 62, 62, 63, 64, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 70, 71, 72, 72, 72, 73, 73, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 75, 75, 75, 75, 75, 74, 74, 73, 73, 73, 72, 71, 71, 70, 70, 69, 68, 68, 67, 67, 66, 66, 65, 65, 65, 65, 65, 64, 65, 65, 65, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, 70, 71, 71, 72, 72, 73, 73, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 74, 74, 74, 74, 73, 73, 73, 72, 72, 72, 71, 71, 70, 70, 70, 69, 69, 69, 68, 68, 68, 68, 68, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, 68, 69, 69, 70, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 86, 87, 88, 89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 100, 101, 102, 102, 102, 103, 103, 103, 103, 103, 102, 102, 102, 101, 101, 100, 100, 99, 99, 98, 97, 96, 96, 95, 94, 94, 93, 92, 92, 91, 90, 90, 89, 89, 88, 88, 87, 87, 87, 86, 86, 85, 85, 85, 84, 83, 83, 82, 81, 81, 80, 79, 78, 77, 75, 74, 73, 72, 70, 69, 68, 66, 65, 63, 62, 61, 59, 58, 57, 55, 54, 53, 52, 52, 51, 50, 50, 49, 49, 49, 49, 49, 49, 50, 50, 51, 51, 52, 53, 54, 55, 56, 57, 58, 59, 59, 60, 61, 62, 62, 63, 63, 64, 64, 64, 64, 64, 63, 63, 62, 61, 60, 59, 58, 56, 55, 53, 52, 50, 48, 46, 45, 43, 41, 39, 37, 35, 33, 31, 29, 28, 26, 24, 22, 21, 19, 18, 16, 15, 14, 13, 12, 11, 10, 10, 9, 8, 8, 8, 8, 7, 7, 8, 8, 8, 8, 9, 9, 10, 11, 11, 12, 13, 14, 15, 16, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29, 30, 30, 31, 31, 32, 32, 32, 32, 32, 31, 31, 30, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 19, 18, 17, 16, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, -0, -0, -1, -1, -2, -2, -3, -3, -4, -5, -6, -7, -7, -8, -9, -10, -11, -13, -14, -15, -16, -17, -18, -19, -20, -21, -21, -22, -23, -23, -24, -24, -25, -25, -25, -25, -25, -24, -24, -24, -23, -23, -22, -22, -21, -20, -19, -18, -17, -17, -16, -15, -14, -13, -12, -11, -10, -10, -9, -8, -8, -7, -7, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -7, -7, -8, -8, -9, -9, -10, -11, -12, -12, -13, -14, -15, -16, -17, -18, -19, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -29, -30, -31, -32, -32, -33, -34, -34, -34, -35, -35, -35, -35, -35, -35, -35, -35, -35, -34, -34, -33, -33, -32, -32, -31, -31, -30, -29, -29, -28, -28, -27, -26, -26, -25, -25, -24, -24, -23, -23, -23, -22, -22, -21, -21, -21, -21, -20, -20, -20, -20, -20, -20, -19, -19, -19, -19, -19, -19, -20, -20, -20, -20, -20, -20, -21, -21, -21, -21, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -22, -21, -21, -21, -21, -20, -20, -20, -20, -20, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -19, -20, -20, -20, -20, -19, -19, -19, -19, -19, -19, -18, -18, -18, -17, -17, -17, -16, -16, -16, -15, -15, -15, -15, -14, -14, -14, -14, -14, -15, -15, -15, -16, -16, -17, -18, -19, -20, -21, -22, -23, -24, -26, -27, -28, -30, -31, -33, -34, -35, -37, -38, -39, -40, -41, -42, -43, -44, -45, -46, -47, -48, -48, -49, -49, -50, -50, -51, -51, -52, -52, -52, -53, -53, -53, -53, -53, -53, -54, -54, -54, -54, -54, -54, -54, -54, -54, -54, -54, -54, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -52, -52, -52, -52, -52, -51, -51, -51, -51, -50, -50, -50, -49, -49, -49, -48, -48, -48, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -47, -48, -48, -49, -49, -50, -50, -51, -51, -52, -53, -53, -54, -54, -55, -55, -56, -57, -57, -57, -58, -58, -59, -59, -59, -60, -60, -60, -60, -61, -61, -61, -61, -61, -62, -62, -62, -62, -63, -63, -63, -64, -64, -65, -65, -66, -66, -67, -67, -68, -69, -69, -70, -70, -71, -71, -72, -72, -73, -73, -73, -73, -73, -73, -73, -73, -73, -72, -72, -71, -70, -69, -68, -67, -66, -65, -64, -62, -61, -59, -57, -56, -54, -52, -51, -49, -47, -46, -44, -42, -40, -39, -37, -36, -34, -33, -32, -30, -29, -28, -27, -26, -25, -24, -24, -23, -23, -22, -22, -22, -22, -22, -22, -23, -23, -24, -24, -25, -26, -27, -28, -29, -31, -32, -34, -35, -37, -39, -41, -43, -45, -47, -49, -51, -54, -56, -58, -61, -63, -66, -68, -71, -73, -75, -77, -80, -82, -84, -85, -87, -89, -90, -91, -92, -93, -94, -95, -95, -95, -96, -95, -95, -95, -95, -94, -93, -93, -92, -91, -90, -89, -89, -88, -87, -86, -85, -84, -84, -83, -82, -82, -81, -81, -80, -80, -79, -79, -79, -78, -78, -78, -78, -77, -77, -77, -76, -76, -75, -75, -74, -74, -73, -72, -72, -71, -70, -69, -68, -67, -67, -66, -65, -64, -63, -63, -62, -62, -61, -61, -61, -61, -61, -61, -61, -62, -62, -63, -64, -65, -66, -67, -68, -70, -71, -72, -74, -75, -76, -78, -79, -80, -81, -83, -84, -84, -85, -86, -86, -87, -87, -87, -87, -87, -87, -87, -87, -86, -86, -85, -84, -83, -82, -82, -81, -80, -79, -78, -77, -75, -74, -74, -73, -72, -71, -70, -69, -69, -68, -68, -68, -68, -68, -68, -68, -69, -69, -70, -71, -72, -73, -74, -75, -76, -78, -79, -81, -82, -84, -86, -87, -89, -90, -92, -93, -95, -96, -97, -98, -99, -100, -101, -102, -102, -102, -103, -103, -103, -103, -103, -102, -102, -101, -101, -100, -99, -98, -97, -96, -95, -94, -93, -92, -91, -89, -88, -87, -86, -85, -84, -83, -82, -81, -80, -79, -78, -77, -76, -76, -75, -75, -74, -74, -74, -74, -74, -74, -74, -74, -74, -74, -75, -75, -76, -76, -77, -78, -78, -79, -80, -81, -82, -83, -84, -85, -86, -87, -88, -88, -89, -90, -91, -92, -93, -94, -94, -95, -96, -96, -97, -97, -97, -98, -98, -98, -98, -98, -98, -98, -98, -98, -98, -98, -98, -97, -97, -97, -97, -96, -96, -95, -95, -94, -94, -93, -93, -92, -91, -91, -90, -90, -89, -88, -88, -87, -87, -86, -86, -85, -85, -84, -84, -84, -83, -83, -83, -83, -83, -83, -83, -83, -83, -84, -84, -84, -85, -85, -85, -86, -86, -87, -87, -87, -88, -88, -89, -89, -90, -90, -91, -91, -92, -92, -93, -93, -94, -94, -95, -95, -96, -97, -97, -98, -99, -99, -100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111, -112, -113, -114, -115, -116, -117, -118, -119, -120, -121, -122, -122, -123, -121, -122, -123, -124, -125, -126, -127, -128, -129, -129, -130, -131, -131, -132, -132, -133, -133, -133, -134, -134, -135, -135, -136, -136, -137, -137, -138, -138, -139, -140, -141, -141, -142, -143, -144, -145, -145, -146, -147, -148, -149, -150, -151, -152, -153, -154, -155, -156, -157, -157, -158, -159, -160, -160, -161, -162, -162, -163, -163, -163, -164, -164, -164, -164, -164, -164, -164, -164, -164, -164, -163, -163, -163, -163, -162, -162, -162, -161, -161, -160, -160, -159, -159, -158, -157, -157, -156, -155, -155, -154, -153, -152, -151, -150, -149, -148, -146, -145, -144, -143, -141, -140, -138, -137, -135, -134, -132, -131, -129, -127, -126, -124, -123, -121, -119, -118, -116, -115, -113, -112, -110, -109, -108, -107, -106, -104, -103, -102, -102, -101, -100, -100, -99, -99, -99, -98, -98, -98, -99, -99, -99, -100, -100, -101, -102, -102, -103, -104, -105, -106, -107, -108, -109, -111, -112, -113, -114, -116, -117, -118, -119, -121, -122, -123, -124, -126, -127, -128, -129, -130, -131, -132, -133, -134, -135, -136, -136, -137, -137, -138, -138, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -139, -139, -140, -140, -141, -142, -143, -144, -145, -146, -147, -148, -150, -151, -153, -154, -156, -158, -159, -161, -163, -164, -166, -168, -169, -171, -172, -173, -175, -176, -177, -178, -179, -180, -181, -181, -182, -182, -182, -183, -183, -183, -183, -183, -182, -182, -182, -181, -181, -180, -180, -179, -179, -178, -178, -177, -176, -176, -175, -174, -174, -173, -173, -172, -172, -171, -171, -171, -170, -170, -170, -170, -170, -170, -170, -170, -170, -170, -171, -171, -171, -172, -173, -173, -174, -175, -176, -177, -178, -179, -180, -181, -182, -183, -185, -186, -187, -188, -190, -191, -193, -194, -195, -197, -198, -199, -200, -202, -203, -204, -205, -206, -207, -208, -209, -210, -211, -211, -212, -212, -213, -213, -213, -213, -213, -213, -213, -213, -213, -212, -212, -211, -210, -210, -209, -208, -207, -206, -205, -204, -202, -201, -200, -199, -197, -196, -194, -193, -191, -190, -188, -187, -185, -184, -182, -181, -179, -178, -176, -175, -174, -172, -171, -170, -169, -168, -167, -166, -165, -164, -163, -163, -162, -161, -161, -160, -160, -159, -159, -159, -158, -158, -158, -158, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -158, -158, -158, -158, -158, -158, -158, -158, -158, -158, -159, -159, -159, -159, -159, -159, -160, -160, -160, -160, -160, -161, -161, -161, -162, -162, -162, -163, -163, -163, -164, -164, -165, -165, -165, -166, -166, -167, -167, -168, -169, -169, -170, -170, -171, -171, -172, -173, -173, -174, -175, -175, -176, -176, -177, -178, -178, -179, -179, -180, -180, -181, -181, -182, -182, -183, -183, -183, -183, -184, -184, -184, -184, -184, -184, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -185, -184, -184, -184, -184, -184, -184, -183, -183, -183, -182, -182, -182, -181, -181, -180, -180, -179, -179, -179, -178, -178, -177, -177, -176, -176, -175, -175, -174, -174, -173, -173, -173, -172, -172, -172, -171, -171, -171, -171, -170, -170, -170, -170, -170, -170, -169, -169, -169, -169, -169, -169, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -168, -169, -169, -169, -170, -170, -171, -171, -172, -173, -174, -174, -175, -176, -177, -177, -178, -179, -180, -180, -181, -182, -182, -183, -183, -184, -185, -185, -185, -186, -186, -187, -187, -187, -188, -188, -188, -188, -189, -189, -189, -189, -189, -190, -190, -190, -190, -190, -190, -190, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -191, -192, -192, -192, -192, -193, -193, -193, -194, -194, -195, -195, -196, -197, -198, -198, -199, -200, -201, -202, -203, -204, -205, -206, -207, -208, -209, -210, -211, -212, -213, -214, -215, -215, -216, -217, -218, -218, -219, -219, -220, -220, -220, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -221, -220, -220, -220, -220, -220, -220, -220, -220, -220, -220, -220, -221, -221, -221, -221, -221, -222, -222, -222, -222, -223, -223, -223, -224, -224, -224, -224, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -224, -224, -224, -224, -223, -223, -223, -222, -222, -221, -221, -221, -220, -220, -220, -219, -219, -219, -218, -218, -218, -217, -217, -217, -217, -217, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -216, -217, -217, -217, -217, -218, -218, -218, -219, -219, -220, -220, -220, -221, -221, -222, -222, -222, -223, -223, -223, -224, -224, -224, -224, -224, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -225, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -224, -225, -225, -225, -225, -226, -226, -227, -228, -228, -229, -230, -231, -232, -233, -234, -235, -236, -237, -238, -239, -240, -241, -243, -244, -245, -246, -247, -248, -248, -249, -250, -251, -251, -252, -252, -253, -253, -253, -253, -253, -253, -253, -253, -253, -252, -252, -251, -251, -250, -249, -248, -248, -247, -246, -245, -244, -243, -241, -240, -239, -238, -237, -236, -235, -234, -233, -232, -232, -231, -230, -230, -229, -229, -229, -228, -228, -228, -228, -228, -229, -229, -229, -230, -230, -230, -231, -231, -232, -232, -233, -234, -234, -235, -235, -236, -236, -237, -237, -237, -238, -238, -238, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -239, -240, -240, -240, -240, -240, -240, -240, -240, -241, -241, -241, -241, -241, -240, -240, -240, -240, -240, -239, -239, -239, -238, -238, -237, -237, -236, -236, -235, -235, -235, -234, -234, -233, -233, -233, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -232, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -231, -232, -232, -232, -233, -233, -233, -234, -234, -234, -235, -235, -235, -236, -236, -236, -236, -237, -237, -237, -237, -237, -237, -237, -236, -236, -236, -236, -236, -235, -235, -235, -235, -234, -234, -234, -234, -234, -234, -234, -234, -234, -235, -235, -236, -236, -237, -237, -238, -239, -240, -241, -242, -243, -244, -245, -246, -247, -248, -249, -250, -251, -252, -253, -253, -254, -254, -255, -255, -255, -256, -256, -256, -255, -255, -255, -255, -254, -254, -253, -253, -252, -251, -250, -250, -249, -248, -247, -247, -246, -245, -244, -244, -243, -243, -242, -242, -241, -241, -241, -240, -240, -240, -240, -240, -240, -241, -241, -241, -242, -242, -243, -243, -244, -244, -245, -246, -246, -247, -247, -248, -249, -249, -250, -250, -251, -251, -252, -252, -252, -253, -253, -253, -253, -253, -253, -254, -254, -254, -254, -254, -254, -254, -254, -254, -254, -254, -254, -254, -255, -255, -255, -255, -256, -256, -256, -257, -257, -258, -258, -259, -259, -259, -260, -260, -261, -261, -262, -262, -262, -263, -263, -263, -264, -264, -264, -264, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -265, -264, -264, -264, -263, -263, -262, -262, -261, -261, -260, -259, -258, -258, -257, -256, -255, -254, -253, -252, -251, -250, -249, -248, -247, -246, -245, -244, -243, -242, -241, -240, -239, -238, -237, -236, -236, -235, -234, -234, -233, -233, -232, -232, -231, -231, -230, -230, -230, -229, -229, -229, -229, -228, -228, -228, -228, -227, -227, -227, -226, -226, -226, -225, -225, -224, -223, -223, -222, -221, -221, -220, -219, -218, -214, -213, -212, -210, -207, -206, -205, -201, -199, -197, -196, -194, -192, -191, -180, -171, -154};
int defaultpulse[8000] = {-97, -97, -97, -97, -97, -97, -98, -98, -98, -98, -98, -98, -98, -99, -99, -99, -99, -100, -100, -100, -100, -101, -101, -101, -101, -102, -102, -102, -103, -103, -103, -103, -104, -104, -104, -105, -105, -105, -106, -106, -106, -107, -107, -107, -107, -108, -108, -108, -109, -109, -109, -109, -110, -110, -110, -110, -111, -111, -111, -111, -111, -112, -112, -112, -112, -112, -112, -112, -113, -113, -113, -113, -113, -113, -113, -113, -113, -114, -114, -114, -114, -114, -114, -114, -114, -115, -115, -115, -115, -115, -115, -116, -116, -116, -116, -116, -117, -117, -117, -117, -118, -118, -118, -118, -119, -119, -119, -119, -120, -120, -120, -121, -121, -121, -122, -122, -122, -122, -123, -123, -123, -124, -124, -124, -125, -125, -125, -125, -126, -126, -126, -126, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -127, -126, -126, -126, -126, -125, -125, -125, -124, -124, -124, -123, -123, -122, -122, -122, -121, -121, -120, -120, -119, -119, -118, -118, -117, -117, -116, -116, -115, -115, -114, -114, -114, -113, -113, -113, -112, -112, -112, -112, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -112, -113, -113, -113, -113, -113, -113, -113, -113, -113, -113, -114, -114, -114, -114, -114, -115, -115, -115, -115, -116, -116, -116, -116, -117, -117, -117, -118, -118, -118, -119, -119, -120, -120, -120, -121, -121, -121, -122, -122, -122, -123, -123, -123, -123, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -123, -123, -123, -122, -122, -121, -121, -120, -120, -119, -119, -118, -118, -117, -116, -116, -115, -114, -114, -113, -113, -112, -112, -111, -111, -111, -110, -110, -110, -109, -109, -109, -109, -109, -109, -109, -109, -109, -109, -110, -110, -110, -110, -111, -111, -111, -112, -112, -112, -113, -113, -113, -114, -114, -114, -115, -115, -115, -116, -116, -116, -116, -117, -117, -117, -117, -118, -118, -118, -118, -118, -119, -119, -119, -119, -119, -119, -120, -120, -120, -120, -120, -121, -121, -121, -121, -121, -121, -122, -122, -122, -122, -122, -122, -123, -123, -123, -123, -123, -123, -123, -124, -124, -124, -124, -124, -124, -124, -124, -125, -125, -125, -125, -125, -125, -125, -125, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -126, -127, -127, -127, -127, -127, -127, -128, -128, -128, -129, -129, -129, -130, -130, -130, -131, -131, -132, -132, -132, -133, -133, -133, -134, -134, -134, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -133, -133, -133, -133, -132, -132, -132, -132, -132, -131, -131, -131, -131, -130, -130, -130, -130, -130, -129, -129, -129, -129, -129, -128, -128, -128, -128, -127, -127, -127, -127, -126, -126, -126, -126, -126, -125, -125, -125, -125, -125, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -125, -125, -125, -125, -126, -126, -126, -126, -127, -127, -127, -128, -128, -129, -129, -129, -130, -130, -130, -131, -131, -131, -132, -132, -132, -132, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -132, -132, -132, -132, -132, -132, -131, -131, -131, -131, -131, -131, -130, -130, -130, -130, -130, -130, -130, -130, -130, -129, -129, -129, -129, -129, -129, -129, -129, -129, -129, -129, -128, -128, -128, -128, -128, -128, -128, -128, -127, -127, -127, -127, -127, -126, -126, -126, -126, -126, -126, -125, -125, -125, -125, -125, -125, -125, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -124, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -126, -126, -126, -126, -126, -126, -126, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -126, -126, -126, -126, -126, -126, -127, -127, -127, -127, -127, -127, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -127, -127, -127, -127, -127, -127, -127, -127, -126, -126, -126, -126, -126, -126, -125, -125, -125, -125, -124, -124, -124, -123, -123, -123, -122, -122, -122, -121, -121, -121, -120, -120, -120, -119, -119, -118, -118, -118, -117, -117, -117, -117, -116, -116, -116, -115, -115, -115, -115, -114, -114, -114, -114, -114, -113, -113, -113, -113, -112, -112, -112, -112, -112, -112, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -111, -112, -112, -112, -112, -113, -113, -113, -114, -114, -114, -115, -115, -116, -116, -117, -117, -118, -118, -119, -119, -120, -120, -121, -121, -122, -123, -123, -124, -124, -125, -125, -126, -126, -127, -127, -128, -128, -129, -129, -129, -130, -130, -130, -131, -131, -131, -131, -131, -132, -132, -132, -132, -132, -132, -132, -132, -132, -131, -131, -131, -131, -131, -131, -131, -130, -130, -130, -130, -129, -129, -129, -129, -129, -128, -128, -128, -128, -128, -127, -127, -127, -127, -127, -127, -126, -126, -126, -126, -126, -126, -126, -126, -126, -127, -127, -127, -127, -127, -127, -128, -128, -128, -129, -129, -129, -130, -130, -130, -131, -131, -132, -132, -132, -133, -133, -133, -134, -134, -134, -135, -135, -135, -135, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -133, -133, -133, -133, -133, -133, -133, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -132, -132, -132, -132, -132, -132, -133, -133, -133, -133, -134, -134, -134, -135, -135, -135, -136, -136, -136, -137, -137, -137, -138, -138, -138, -138, -139, -139, -139, -139, -139, -139, -140, -140, -140, -140, -140, -140, -140, -140, -140, -140, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -138, -138, -138, -138, -138, -138, -137, -137, -137, -137, -137, -137, -136, -136, -136, -136, -136, -136, -136, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -136, -136, -136, -136, -136, -136, -136, -136, -137, -137, -137, -137, -137, -137, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -137, -137, -137, -137, -136, -136, -136, -136, -135, -135, -135, -134, -134, -134, -134, -133, -133, -133, -133, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -133, -133, -133, -133, -133, -133, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -136, -136, -136, -136, -136, -136, -136, -137, -137, -137, -137, -137, -137, -137, -137, -137, -137, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -137, -137, -137, -136, -136, -136, -135, -135, -135, -134, -134, -134, -133, -133, -133, -132, -132, -132, -132, -131, -131, -131, -131, -131, -131, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -131, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -133, -133, -133, -133, -133, -133, -133, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -134, -134, -134, -133, -133, -133, -133, -132, -132, -132, -132, -131, -131, -131, -131, -130, -130, -130, -130, -129, -129, -129, -129, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -128, -129, -129, -129, -130, -130, -130, -131, -131, -132, -132, -133, -133, -134, -134, -135, -135, -136, -136, -137, -137, -138, -138, -139, -139, -140, -140, -140, -141, -141, -141, -141, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -140, -140, -140, -140, -140, -140, -140, -140, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -137, -137, -137, -137, -137, -137, -137, -137, -136, -136, -136, -136, -136, -135, -135, -135, -135, -134, -134, -134, -134, -133, -133, -133, -133, -133, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -132, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -136, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -134, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -134, -134, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -133, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -136, -136, -136, -136, -137, -137, -137, -137, -137, -137, -137, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -137, -137, -137, -137, -137, -137, -137, -137, -137, -137, -136, -136, -136, -136, -136, -136, -136, -135, -135, -135, -135, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -135, -136, -136, -136, -136, -136, -136, -136, -136, -137, -137, -137, -137, -137, -137, -137, -137, -137, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -137, -137, -137, -137, -137, -137, -137, -138, -138, -138, -138, -138, -138, -138, -139, -139, -139, -140, -140, -140, -141, -141, -142, -142, -142, -143, -143, -144, -144, -145, -146, -146, -147, -147, -147, -148, -148, -149, -149, -150, -150, -150, -151, -151, -151, -152, -152, -152, -152, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -153, -152, -152, -152, -152, -151, -151, -151, -150, -150, -150, -149, -149, -148, -148, -147, -147, -146, -146, -145, -144, -144, -143, -143, -142, -141, -141, -140, -140, -139, -139, -138, -138, -137, -137, -136, -136, -135, -135, -135, -134, -134, -134, -133, -133, -133, -132, -132, -132, -131, -131, -131, -131, -130, -130, -130, -129, -129, -129, -128, -128, -128, -128, -127, -127, -127, -126, -126, -126, -126, -126, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -125, -126, -126, -126, -126, -127, -127, -127, -128, -128, -128, -129, -129, -129, -130, -130, -131, -131, -131, -132, -132, -132, -133, -133, -133, -134, -134, -134, -135, -135, -135, -136, -136, -136, -136, -137, -137, -137, -137, -137, -138, -138, -138, -138, -139, -139, -139, -139, -140, -140, -140, -140, -140, -141, -141, -141, -141, -141, -141, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -141, -141, -141, -141, -141, -140, -140, -140, -140, -139, -139, -139, -139, -138, -138, -138, -137, -137, -136, -136, -136, -135, -135, -135, -134, -134, -134, -133, -133, -133, -132, -132, -132, -131, -131, -131, -131, -131, -130, -130, -130, -130, -130, -130, -130, -130, -130, -131, -131, -131, -131, -132, -132, -132, -133, -133, -133, -134, -134, -135, -135, -136, -136, -136, -137, -137, -138, -138, -139, -139, -140, -140, -140, -141, -141, -141, -142, -142, -142, -142, -143, -143, -143, -143, -143, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -145, -145, -145, -145, -145, -146, -146, -146, -146, -146, -147, -147, -147, -147, -148, -148, -148, -148, -149, -149, -149, -149, -149, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -149, -149, -149, -149, -149, -149, -148, -148, -148, -148, -147, -147, -147, -147, -147, -146, -146, -146, -146, -145, -145, -145, -145, -145, -144, -144, -144, -144, -144, -143, -143, -143, -143, -143, -143, -142, -142, -142, -142, -142, -142, -142, -142, -142, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -141, -140, -140, -140, -140, -140, -140, -140, -140, -140, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -139, -139, -139, -139, -139, -139, -139, -139, -139, -140, -140, -140, -140, -141, -141, -141, -141, -142, -142, -142, -142, -143, -143, -143, -143, -144, -144, -144, -144, -144, -145, -145, -145, -145, -145, -145, -145, -145, -146, -146, -146, -146, -146, -146, -146, -146, -146, -146, -146, -146, -146, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -146, -146, -146, -146, -146, -146, -145, -145, -145, -145, -145, -145, -145, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -145, -145, -145, -145, -145, -145, -146, -146, -146, -146, -147, -147, -147, -147, -148, -148, -148, -148, -148, -148, -149, -149, -149, -149, -149, -149, -149, -149, -149, -148, -148, -148, -148, -148, -147, -147, -147, -146, -146, -146, -145, -145, -144, -144, -143, -143, -142, -142, -141, -141, -140, -140, -139, -139, -138, -138, -138, -137, -137, -136, -136, -136, -135, -135, -135, -135, -134, -134, -134, -134, -134, -134, -134, -135, -135, -135, -135, -136, -136, -136, -137, -137, -138, -138, -138, -139, -140, -140, -141, -141, -142, -142, -143, -143, -144, -144, -144, -145, -145, -146, -146, -146, -146, -147, -147, -147, -147, -147, -147, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -149, -149, -149, -149, -149, -149, -149, -149, -149, -149, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -151, -151, -151, -151, -151, -151, -151, -151, -152, -152, -152, -152, -152, -152, -153, -153, -153, -153, -153, -153, -153, -153, -153, -154, -154, -154, -154, -154, -154, -154, -154, -154, -154, -153, -153, -153, -153, -153, -153, -152, -152, -152, -152, -151, -151, -151, -151, -150, -150, -150, -149, -149, -149, -149, -148, -148, -148, -148, -148, -148, -148, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -146, -146, -146, -146, -146, -146, -146, -145, -145, -145, -145, -145, -145, -145, -144, -144, -144, -144, -144, -144, -144, -144, -144, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -144, -144, -144, -144, -144, -144, -145, -145, -145, -145, -145, -146, -146, -146, -146, -147, -147, -147, -148, -148, -148, -148, -149, -149, -149, -149, -149, -149, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -149, -149, -149, -149, -149, -148, -148, -148, -147, -147, -147, -146, -146, -146, -145, -145, -145, -144, -144, -143, -143, -143, -142, -142, -142, -141, -141, -141, -141, -140, -140, -140, -140, -140, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -138, -139, -139, -139, -139, -139, -140, -140, -140, -140, -141, -141, -141, -142, -142, -142, -143, -143, -143, -144, -144, -144, -145, -145, -145, -146, -146, -146, -147, -147, -148, -148, -148, -149, -149, -149, -150, -150, -151, -151, -151, -152, -152, -152, -153, -153, -153, -154, -154, -154, -155, -155, -155, -156, -156, -156, -156, -156, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -157, -156, -156, -156, -156, -155, -155, -155, -155, -154, -154, -153, -153, -153, -152, -152, -151, -151, -150, -150, -149, -148, -148, -147, -147, -146, -145, -145, -144, -144, -143, -143, -142, -142, -141, -141, -141, -140, -140, -140, -140, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -139, -140, -140, -140, -140, -140, -141, -141, -141, -141, -142, -142, -142, -142, -143, -143, -143, -143, -144, -144, -144, -144, -145, -145, -145, -145, -146, -146, -146, -147, -147, -147, -147, -148, -148, -148, -149, -149, -149, -149, -150, -150, -150, -150, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -150, -150, -150, -149, -149, -149, -149, -148, -148, -148, -147, -147, -147, -146, -146, -146, -145, -145, -145, -145, -144, -144, -144, -144, -144, -144, -144, -143, -143, -143, -143, -143, -143, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -145, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -144, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -143, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -142, -143, -143, -143, -143, -144, -144, -144, -144, -145, -145, -145, -146, -146, -146, -146, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -146, -146, -146, -146, -145, -145, -145, -144, -144, -144, -143, -143, -142, -142, -142, -142, -141, -141, -141, -141, -141, -140, -140, -140, -140, -140, -140, -140, -141, -141, -141, -141, -141, -141, -142, -142, -142, -143, -143, -143, -143, -144, -144, -144, -144, -145, -145, -145, -145, -146, -146, -146, -146, -146, -147, -147, -147, -147, -147, -147, -147, -147, -147, -147, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -148, -149, -149, -149, -149, -149, -149, -149, -149, -149, -149, -149, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -150, -151, -151, -151, -151, -151, -151, -151, -151, -151, -151, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -152, -151, -151, -151, -151, -151, -150, -150, -150, -149, -149, -148, -147, -147, -146, -145, -144, -143, -141, -140, -138, -137, -135, -133, -130, -128, -125, -123, -120, -117, -113, -110, -106, -102, -98, -94, -90, -85, -81, -76, -71, -66, -60, -55, -50, -44, -39, -33, -27, -22, -16, -11, -5, -0, 5, 10, 15, 20, 25, 29, 33, 37, 41, 44, 48, 51, 53, 56, 58, 60, 61, 62, 63, 64, 64, 65, 64, 64, 63, 62, 61, 60, 58, 57, 55, 53, 50, 48, 45, 43, 40, 37, 34, 31, 28, 25, 22, 19, 16, 13, 9, 6, 3, -0, -3, -6, -9, -12, -15, -18, -21, -24, -27, -30, -32, -35, -37, -40, -42, -45, -47, -49, -52, -54, -56, -58, -60, -62, -64, -66, -68, -69, -71, -73, -75, -76, -78, -79, -81, -82, -83, -85, -86, -87, -88, -90, -91, -92, -93, -94, -95, -96, -97, -98, -99, -100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -110, -111, -112, -113, -114, -115, -116, -117, -118, -119, -120, -120, -121, -122, -123, -124, -125, -126, -126, -127, -128, -129, -129, -130, -131, -131, -132, -133, -133, -134, -134, -134, -135, -135, -135, -136, -136, -136, -136, -136, -137, -137, -137, -137, -137, -136, -136, -136, -136, -136, -136, -135, -135, -135, -134, -134, -134, -133, -132, -132, -131, -131, -130, -129, -128, -127, -126, -125, -124, -123, -122, -120, -119, -118, -116, -115, -113, -111, -109, -108, -106, -104, -101, -99, -97, -94, -92, -89, -86, -83, -80, -76, -72, -68, -64, -59, -54, -49, -44, -38, -31, -24, -17, -9, -1, 8, 17, 27, 37, 48, 60, 72, 84, 97, 111, 125, 139, 154, 170, 186, 202, 218, 235, 252, 270, 287, 305, 322, 340, 358, 375, 392, 409, 426, 442, 458, 473, 488, 502, 515, 528, 539, 550, 560, 569, 576, 583, 589, 593, 596, 599, 600, 600, 598, 596, 592, 588, 582, 575, 568, 559, 550, 539, 528, 517, 505, 492, 479, 465, 451, 437, 423, 409, 395, 381, 367, 353, 339, 326, 313, 301, 289, 277, 266, 256, 246, 236, 227, 219, 211, 204, 198, 192, 186, 181, 176, 172, 168, 164, 161, 158, 155, 153, 150, 148, 146, 144, 142, 139, 137, 135, 133, 130, 128, 125, 123, 120, 117, 114, 110, 107, 103, 99, 95, 91, 87, 83, 79, 74, 70, 65, 61, 56, 52, 47, 43, 38, 34, 30, 26, 22, 18, 14, 11, 7, 4, 1, -1, -4, -6, -9, -11, -12, -14, -15, -16, -17, -17, -18, -18, -18, -18, -17, -16, -16, -15, -13, -12, -10, -9, -7, -5, -3, -1, 1, 4, 6, 8, 10, 13, 15, 17, 19, 21, 23, 25, 27, 28, 30, 31, 33, 34, 35, 35, 36, 37, 37, 37, 37, 37, 37, 37, 36, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 16, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 16, 17, 17, 18, 19, 20, 21, 22, 24, 25, 26, 28, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42, 44, 46, 47, 49, 52, 54, 56, 59, 62, 66, 70, 75, 80, 86, 92, 100, 108, 118, 128, 141, 154, 169, 186, 205, 226, 249, 275, 303, 333, 366, 403, 442, 484, 530, 579, 631, 687, 746, 809, 875, 945, 1018, 1095, 1174, 1257, 1343, 1431, 1522, 1615, 1710, 1806, 1904, 2003, 2103, 2203, 2303, 2402, 2500, 2597, 2693, 2786, 2878, 2966, 3051, 3133, 3211, 3285, 3355, 3421, 3481, 3537, 3588, 3635, 3676, 3711, 3742, 3768, 3788, 3804, 3815, 3821, 3823, 3820, 3813, 3802, 3787, 3769, 3747, 3722, 3694, 3664, 3631, 3596, 3559, 3521, 3480, 3439, 3396, 3353, 3309, 3264, 3219, 3174, 3128, 3083, 3038, 2993, 2948, 2904, 2861, 2818, 2775, 2733, 2692, 2651, 2611, 2572, 2533, 2495, 2457, 2419, 2382, 2345, 2308, 2272, 2236, 2200, 2164, 2128, 2092, 2056, 2020, 1984, 1948, 1912, 1876, 1839, 1803, 1767, 1730, 1694, 1658, 1622, 1586, 1550, 1514, 1479, 1443, 1409, 1374, 1340, 1307, 1273, 1241, 1209, 1178, 1147, 1117, 1088, 1059, 1031, 1004, 978, 953, 928, 905, 882, 860, 839, 819, 800, 781, 763, 746, 730, 715, 700, 686, 672, 659, 647, 635, 624, 613, 603, 593, 583, 573, 564, 555, 547, 538, 530, 522, 514, 506, 498, 491, 483, 476, 468, 461, 454, 447, 441, 434, 427, 421, 414, 408, 402, 396, 391, 385, 380, 374, 369, 364, 359, 355, 350, 345, 341, 336, 332, 328, 324, 319, 315, 311, 307, 303, 299, 294, 290, 286, 281, 277, 272, 268, 263, 258, 254, 249, 244, 239, 234, 229, 224, 219, 214, 209, 204, 199, 195, 190, 185, 181, 177, 173, 169, 165, 161, 157, 154, 151, 148, 145, 143, 140, 138, 136, 134, 132, 130, 129, 128, 126, 125, 124, 123, 123, 122, 121, 120, 120, 119, 119, 118, 117, 117, 116, 115, 115, 114, 113, 112, 111, 110, 109, 108, 107, 106, 105, 104, 103, 101, 100, 99, 98, 96, 95, 94, 92, 91, 90, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 79, 78, 77, 77, 77, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 78, 78, 79, 79, 79, 80, 80, 80, 81, 81, 82, 82, 82, 83, 83, 84, 84, 84, 84, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 85, 85, 85, 85, 84, 84, 84, 83, 83, 83, 82, 82, 81, 81, 80, 80, 79, 78, 78, 77, 77, 76, 76, 75, 75, 74, 74, 73, 73, 73, 72, 72, 72, 72, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 73, 73, 73, 74, 74, 74, 75, 75, 75, 75, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 76, 76, 76, 76, 75, 75, 75, 74, 74, 74, 73, 73, 72, 72, 71, 71, 71, 70, 70, 70, 70, 69, 69, 69, 69, 69, 69, 69, 69, 70, 70, 70, 70, 71, 71, 72, 72, 72, 73, 73, 74, 74, 75, 75, 76, 76, 77, 77, 78, 78, 78, 79, 79, 79, 79, 80, 80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 84, 84, 84, 85, 85, 85, 85, 86, 86, 86, 86, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 86, 86, 86, 86, 85, 85, 85, 84, 84, 83, 83, 82, 82, 81, 81, 80, 80, 79, 79, 78, 78, 78, 77, 77, 77, 76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 78, 78, 79, 80, 80, 81, 82, 83, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 94, 95, 96, 96, 97, 97, 98, 98, 99, 99, 99, 99, 99, 99, 99, 99, 99, 98, 98, 98, 97, 97, 97, 96, 96, 96, 95, 95, 94, 94, 94, 93, 93, 93, 93, 92, 92, 92, 92, 92, 92, 91, 91, 91, 91, 91, 91, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 92, 92, 92, 91, 91, 91, 90, 90, 89, 89, 88, 88, 87, 87, 87, 86, 86, 85, 85, 84, 84, 83, 83, 83, 82, 82, 82, 81, 81, 81, 81, 80, 80, 80, 80, 80, 80, 80, 80, 81, 81, 81, 82, 82, 82, 83, 84, 84, 85, 86, 86, 87, 88, 89, 90, 91, 91, 92, 93, 94, 95, 95, 96, 97, 97, 98, 99, 99, 100, 100, 100, 101, 101, 101, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 101, 101, 101, 101, 101, 101, 100, 100, 100, 100, 99, 99, 99, 99, 99, 98, 98, 98, 98, 97, 97, 97, 96, 96, 96, 96, 95, 95, 95, 94, 94, 94, 94, 94, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 97, 97, 97, 97, 97, 97, 97, 97, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, 100, 100, 100, 100, 101, 101, 102, 102, 102, 103, 103, 103, 104, 104, 105, 105, 105, 106, 106, 106, 107, 107, 107, 108, 108, 108, 108, 108, 109, 109, 109, 109, 109, 109, 108, 108, 108, 108, 108, 107, 107, 107, 106, 106, 105, 105, 105, 104, 104, 103, 103, 103, 102, 102, 102, 101, 101, 101, 101, 100, 100, 100, 100, 100, 99, 99, 99, 99, 99, 99, 99, 99, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 99, 99, 99, 99, 99, 100, 100, 100, 101, 101, 101, 102, 102, 103, 103, 104, 104, 105, 105, 106, 106, 107, 107, 108, 108, 109, 109, 109, 110, 110, 110, 110, 111, 111, 111, 111, 111, 111, 111, 111, 111, 110, 110, 110, 110, 109, 109, 109, 109, 108, 108, 108, 108, 107, 107, 107, 107, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 105, 105, 105, 105, 105, 105, 105, 104, 104, 104, 104, 103, 103, 103, 102, 102, 101, 101, 101, 100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 96, 95, 95, 94, 94, 93, 93, 93, 92, 92, 92, 91, 91, 91, 91, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 92, 92, 92, 93, 93, 94, 94, 95, 95, 96, 97, 97, 98, 98, 99, 99, 100, 101, 101, 102, 102, 102, 103, 103, 103, 104, 104, 104, 104, 104, 104, 104, 104, 104, 103, 103, 103, 103, 102, 102, 101, 101, 100, 100, 99, 99, 98, 97, 97, 96, 96, 95, 95, 95, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 95, 95, 95, 96, 96, 97, 97, 98, 99, 99, 100, 100, 101, 101, 102, 102, 102, 103, 103, 103, 104, 104, 104, 104, 104, 105, 105, 105, 105, 105, 105, 105, 104, 104, 104, 104, 104, 104, 103, 103, 103, 103, 102, 102, 102, 102, 101, 101, 101, 100, 100, 100, 99, 99, 99, 98, 98, 98, 97, 97, 96, 96, 95, 95, 95, 94, 94, 93, 93, 92, 91, 91, 91, 90, 90, 89, 89, 88, 88, 88, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, 88, 88, 88, 89, 89, 89, 90, 90, 90, 91, 91, 91, 92, 92, 92, 93, 93, 94, 94, 94, 95, 95, 96, 96, 97, 97, 98, 98, 98, 99, 99, 100, 100, 101, 101, 101, 102, 102, 102, 102, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 102, 102, 102, 101, 101, 101, 100, 99, 99, 98, 98, 97, 96, 95, 95, 94, 93, 92, 92, 91, 90, 89, 89, 88, 87, 87, 86, 85, 85, 84, 84, 83, 83, 82, 82, 81, 81, 80, 80, 79, 79, 79, 78, 78, 78, 78, 78, 77, 77, 77, 77, 77, 77, 77, 78, 78, 78, 78, 79, 79, 80, 80, 81, 82, 82, 83, 84, 85, 86, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 96, 97, 98, 99, 99, 100, 100, 101, 101, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 101, 101, 101, 100, 100, 99, 99, 99, 98, 98, 97, 97, 96, 96, 95, 95, 95, 94, 94, 93, 93, 93, 93, 92, 92, 92, 92, 91, 91, 91, 91, 90, 90, 90, 90, 90, 90, 89, 89, 89, 89, 89, 89, 89, 88, 88, 88, 88, 88, 88, 88, 88, 87, 87, 87, 87, 87, 87, 86, 86, 86, 86, 85, 85, 85, 85, 84, 84, 84, 83, 83, 83, 82, 82, 82, 81, 81, 81, 81, 80, 80, 80, 80, 80, 80, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 78, 78, 78, 78, 78, 77, 77, 77, 77, 76, 76, 76, 75, 75, 75, 75, 75, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 78, 78, 78, 79, 79, 80, 80, 81, 81, 81, 82, 82, 83, 83, 84, 84, 84, 85, 85, 86, 86, 86, 86, 86, 87, 87, 87, 87, 87, 87, 87, 87, 87, 86, 86, 86, 86, 86, 85, 85, 85, 84, 84, 84, 83, 83, 83, 82, 82, 82, 81, 81, 81, 81, 80, 80, 80, 79, 79, 79, 79, 78, 78, 78, 77, 77, 77, 77, 76, 76, 76, 75, 75, 75, 74, 74, 73, 73, 72, 72, 71, 71, 70, 70, 69, 69, 68, 68, 68, 67, 67, 66, 66, 66, 66, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 67, 67, 67, 67, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, 69, 68, 68, 68, 68, 68, 67, 67, 67, 66, 66, 65, 65, 64, 64, 64, 63, 63, 62, 62, 62, 61, 61, 61, 61, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, 61, 61, 62, 62, 62, 62, 63, 63, 63, 64, 64, 64, 64, 65, 65, 65, 65, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, 67, 67, 66, 66, 66, 66, 66, 65, 65, 65, 64, 64, 64, 63, 63, 62, 62, 62, 61, 61, 61, 60, 60, 60, 60, 59, 59, 59, 59, 59, 59, 59, 59, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 57, 57, 57, 56, 56, 56, 55, 55, 54, 54, 53, 52, 52, 51, 51, 50, 49, 49, 48, 47, 47, 46, 46, 45, 45, 44, 44, 43, 43, 42, 42, 42, 41, 41, 41, 41, 41, 41, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 42, 42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 49, 49, 50, 50, 51, 52, 52, 53, 53, 54, 54, 55, 55, 55, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 59, 59, 59, 59, 59, 58, 58, 58, 58, 57, 57, 56, 56, 55, 55, 54, 54, 53, 53, 52, 52, 51, 50, 50, 49, 49, 48, 48, 47, 47, 47, 46, 46, 46, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 46, 46, 46, 46, 45, 45, 45, 44, 44, 44, 44, 44, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 47, 47, 47, 47, 46, 46, 45, 45, 44, 44, 43, 42, 42, 41, 40, 39, 39, 38, 37, 36, 35, 35, 34, 33, 33, 32, 32, 31, 31, 30, 30, 30, 30, 29, 29, 30, 30, 30, 30, 30, 31, 31, 32, 32, 33, 33, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 39, 39, 39, 38, 38, 37, 37, 36, 36, 35, 34, 34, 33, 32, 31, 31, 30, 29, 28, 28, 27, 26, 26, 25, 24, 24, 23, 23, 22, 22, 22, 21, 21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 18, 18, 17, 17, 17, 16, 16, 15, 15, 14, 14, 13, 12, 12, 11, 11, 11, 10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 14, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 18, 18, 18, 18, 18, 17, 17, 17, 17, 16, 16, 15, 15, 15, 14, 14, 13, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6, 5, 5, 4, 3, 3, 2, 2, 1, 1, 0, -0, -1, -1, -2, -2, -3, -3, -3, -4, -4, -4, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -4, -4, -4, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 9, 9, 9, 8, 8, 7, 7, 6, 6, 5, 4, 4, 3, 3, 2, 2, 1, 0, -0, -1, -1, -2, -2, -3, -3, -4, -4, -4, -5, -5, -5, -5, -5, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -8, -8, -8, -8, -8, -8, -8, -9, -9, -9, -9, -9, -10, -10, -10, -10, -11, -11, -11, -12, -12, -12, -12, -13, -13, -13, -13, -13, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -14, -13, -13, -13, -13, -13, -13, -13, -13, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -12, -13, -13, -13, -13, -13, -14, -14, -14, -14, -15, -15, -15, -16, -16, -16, -17, -17, -17, -18, -18, -18, -18, -18, -16, -16, -16, -16, -16, -16, -16, -16, -15, -15, -15, -15, -15, -14, -14, -14, -14, -13, -13, -13, -12, -12, -12, -11, -11, -11, -10, -10, -10, -10, -9, -9, -9, -9, -9, -9, -9, -9, -8, -8, -8, -8, -9, -9, -9, -9, -9, -9, -9, -9, -10, -10, -10, -10, -10, -11, -11, -11, -11, -12, -12, -12, -12, -12, -13, -13, -13, -13, -13, -14, -14, -14, -14, -14, -15, -15, -15, -15, -15, -16, -16, -16, -16, -16, -17, -17, -17, -17, -17, -18, -18, -18, -18, -18, -18, -19, -19, -19, -19, -19, -19, -19, -19, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -21, -21, -21, -21, -21, -21, -21, -22, -22, -22, -22, -23, -23, -23, -24, -24, -24, -25, -25, -25, -26, -26, -27, -27, -27, -28, -28, -28, -29, -29, -29, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -29, -29, -29, -29, -29, -29, -29, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -28, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -30, -30, -30, -30, -31, -31, -31, -31, -32, -32, -32, -33, -33, -33, -34, -34, -34, -35, -35, -35, -36, -36, -36, -37, -37, -37, -38, -38, -38, -38, -38, -39, -39, -39, -39, -40, -40, -40, -40, -40, -40, -41, -41, -41, -41, -41, -41, -41, -41, -41, -41, -41, -42, -42, -42, -42, -42, -42, -42, -42, -41, -41, -41, -41, -41, -41, -41, -41, -41, -40, -40, -40, -40, -40, -39, -39, -39, -39, -38, -38, -38, -38, -37, -37, -37, -36, -36, -36, -35, -35, -35, -34, -34, -34, -34, -33, -33, -33, -33, -33, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -32, -33, -33, -33, -33, -33, -33, -34, -34, -34, -34, -34, -35, -35, -35, -35, -35, -36, -36, -36, -36, -36, -36, -36, -37, -37, -37, -37, -37, -38, -38, -38, -38, -38, -39, -39, -39, -40, -40, -40, -41, -41, -41, -42, -42, -43, -43, -44, -44, -44, -45, -45, -46, -46, -46, -47, -47, -47, -48, -48, -48, -48, -48, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -48, -48, -48, -48, -48, -48, -47, -47, -47, -47, -46, -46, -46, -45, -45, -45, -45, -44, -44, -44, -44, -43, -43, -43, -43, -42, -42, -42, -42, -42, -41, -41, -41, -41, -40, -40, -40, -40, -39, -39, -39, -39, -38, -38, -38, -38, -37, -37, -37, -37, -37, -36, -36, -36, -36, -36, -36, -35, -35, -35, -35, -35, -35, -35, -35, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -34, -35, -35, -35, -35, -35, -35, -36, -36, -36, -36, -37, -37, -37, -38, -38, -39, -39, -39, -40, -40, -41, -41, -41, -42, -42, -42, -43, -43, -43, -44, -44, -44, -44, -44, -44, -44, -45, -45, -45, -45, -45, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -43, -43, -43, -43, -43, -43, -43, -43, -43, -43, -43, -43, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -44, -45, -45, -45, -45, -45, -45, -45, -46, -46, -46, -46, -46, -46, -46, -47, -47, -47, -47, -47, -47, -48, -48, -48, -48, -48, -48, -48, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -49, -48, -48, -48, -48, -48, -48, -48, -48, -47, -47, -47, -47, -47, -47, -47, -47, -47, -46, -46, -46, -46, -46, -46, -47, -47, -47, -47, -47, -47, -47, -48, -48, -48, -48, -49, -49, -49, -49, -50, -50, -50, -51, -51, -51, -51, -51, -52, -52, -52, -52, -52, -52, -52, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -53, -54, -54, -54, -54, -54, -54, -54, -54, -54, -55, -55, -55, -55, -55, -55, -55, -55, -55, -55, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -57, -57, -57, -57, -57, -57, -58, -58, -58, -58, -58, -58, -59, -59, -59, -59, -59, -59, -59, -60, -60, -60, -60, -60, -60, -60, -60, -59, -59, -59, -59, -59, -59, -59, -58, -58, -58, -58, -58, -58, -57, -57, -57, -57, -57, -57, -57, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -56, -57, -57, -57, -57, -57, -57, -57, -57, -57, -57, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -58, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -60, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -59, -60, -60, -60, -60, -60, -60, -61, -61, -61, -61, -62, -62, -62, -63, -63, -63, -64, -64, -64, -65, -65, -65, -66, -66, -66, -66, -67, -67, -67, -67, -68, -68, -68, -68, -68, -69, -69, -69, -69, -69, -69, -70, -70, -70, -70, -70, -70, -70, -70, -70, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -71, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -72, -73, -73, -73, -73, -73, -73, -73, -73, -74, -74, -74, -74, -74, -74, -74, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -75, -74, -74, -74, -74, -74, -73, -73, -73, -72, -72, -72, -71, -71, -70, -70, -70, -69, -69, -69, -68, -68, -68, -68, -67, -67, -67, -67, -67, -67, -67, -67, -67, -68, -68, -68, -68, -69, -69, -70, -70, -70, -71, -71, -72, -72, -73, -73, -74, -74, -75, -75, -76, -76, -77, -77, -77, -78, -78, -78, -78, -78, -78, -78, -79, -79, -78, -78, -78, -78, -78, -78, -77, -77, -77, -76, -76, -75, -75, -74, -74, -74, -73, -73, -72, -72, -71, -71, -71, -70, -70, -70, -69, -69, -69, -69, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -68, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -67, -68, -68, -68, -68, -68, -68, -68, -69, -69, -69, -69, -69, -70, -70, -70, -70, -71, -71, -71, -71, -72, -72, -72, -73, -73, -73, -73, -74, -74, -74, -74, -75, -75, -75, -75, -76, -76, -76, -76, -77, -77, -77, -78, -78, -78, -79, -79, -80, -80, -80, -81, -81, -81, -82, -82, -83, -83, -83, -84, -84, -84, -85, -85, -85, -85, -85, -86, -86, -86, -86, -86, -86, -86, -86, -86, -85, -85, -85, -85, -85, -84, -84, -84, -84, -83, -83, -83, -82, -82, -82, -82, -81, -81, -81, -80, -80, -80, -80, -79, -79, -79, -79, -78, -78, -78, -77, -77, -77, -77, -76, -76, -75, -75, -75, -74, -74, -74, -73, -73, -72, -72, -72, -71, -71, -70, -70, -70, -70, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -69, -70, -70, -70, -71, -71, -72, -72, -73, -73, -74, -74, -75, -75, -76, -77, -77, -78, -78, -79, -79, -80, -80, -81, -81, -82, -82, -83, -83, -83, -84, -84, -84, -85, -85, -85, -85, -85, -85, -86, -86, -86, -86, -86, -86, -87, -87, -87, -87, -88, -87, -87, -85, -85, -81, -79, -76, -72, -69, -67, -64, -64, -59, -57, -55, -50, -46, -40, -37, -29, -26, -25, -18, -13, -8, -7, -5, -5, -5, -5, -5, -4, -2, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		fprintf (stderr,"No pulsefile called 'pulseshape' found; using default 8000 shape\nConsider after running with -e 16 : cp shape pulseshape\n");

//                 PEAK=10; eval `cat shape | awk '{print $1,$2/(1+($1>4150)*(($1-4000)/150)^2)/(1+($1<3900)*((-$1+3900)/100)^2),$3/(1+($1>4150)*(($1-4000)/150)^2)/(1+($1<3900)*((-$1+3900)/100)^2)}'  | stat.sh | awk '{printf\"m1=%g;m2=%g;s1=%g;s2=%g\n",$6,$10,$7,$11}'`; cat shape | awk '{print ($1+'$PEAK')%8000,(($2/(1+($1>4150)*(($1-4000)/150)^2)/(1+($1<3900)*((-$1+3900)/100)^2))-('$m1')),(($3/(1+($1>4150)*(($1-4000)/150)^2)/(1+($1<3900)*((-$1+3900)/100)^2)-('$m2')))}'  | sort -g > pulseshape
		for (int j=0; j < NN ; j++) 
{
ps[j] = defaultpulse[j];
pst[j] = svalue?defaulttock[j]:defaultpulse[j];
}
	} 
	else
	{
		for (int p=0; p < NN ; p++) 
		{
			if (fscanf(pulsefile,"%g %g %g", &row,&val,&val2) != 3)
			{
				fprintf (stderr,"pulsefile (pulseshape) should be 1 .0223 0.656 (n tick tock) with %d rows, where peak is at %d\nthe index in column 1 is ignored",NN,NN/2);
				result++;
			}

			ps[p]=0.0;
			ps[p] = val;
			pst[p] = val2;
		}
		fclose(pulsefile);
	}
	return result;
}
