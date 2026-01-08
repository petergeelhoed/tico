/**
 * Read and parse a wave file
 * PG
 * assumes 21600bph and 48000 wav
 **/
#include "mydefs.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// WAVE file header format
struct HEADER
{
    unsigned char riff[4];             // RIFF string
    unsigned int overall_size;         // overall size of file in bytes
    unsigned char wave[4];             // WAVE string
    unsigned char fmt_chunk_marker[4]; // fmt string with trailing null char
    unsigned int length_of_fmt;        // length of the format data
    unsigned int format_type; // format type. 1-PCM, 3- IEEE float, 6 - 8bit A
                              // law, 7 - 8bit mu law
    unsigned int channels;    // no.of channels
    unsigned int sample_rate; // sampling rate (blocks per second)
    unsigned int byterate;    // SampleRate * NumChannels * BitsPerSample/8
    unsigned int block_align; // NumChannels * BitsPerSample/8
    unsigned int bits_per_sample; // bits per sample, 8- 8bits, 16- 16 bits etc
    unsigned char data_chunk_header[4]; // DATA string or FLLR string
    unsigned int data_size; // NumSamples * NumChannels * BitsPerSample/8 - size
                            // of the next chunk that will be read
};

unsigned char buffer4[4];
unsigned char buffer2[2];

struct HEADER header;

struct HEADER readheader();
int openfiles(FILE** tickfile,
              FILE** tockfile,
              FILE** corfile,
              FILE** rawfile,
              FILE* pulsefile,
              FILE** tickavg,
              FILE** teethfile,
              double* ps,
              double* pst,
              int jvalue,
              int wvalue,
              int svalue,
              int NN);

int main()
{
    int rvalue = 0;
    int qvalue = 4000;
    int NN = 8000;
    opterr = 0;
    int read = 0;
    opterr = 0;

    struct HEADER header = readheader();
    // read each sample from data chunk if PCM
    if (header.format_type == 1)
    { // PCM
        long i = 0;
        unsigned char lsb[1];
        signed char msb[1];

        for (i = 44; i <= 100; i++)
        {
            // skip shit
            read = fread(lsb, sizeof(lsb), 1, stdin);
            read += fread(msb, sizeof(msb), 1, stdin);
        }

        long length =
            header.overall_size / 2 + qvalue - rvalue * NN - 100 - NN * 2;

        // loop the entire file
        while (i < length)
        {
            read = fread(lsb, sizeof(lsb), 1, stdin);
            read += fread(msb, sizeof(msb), 1, stdin);
            if (read == 2)
            {
                i++;
                printf("%d\n", (msb[0] << 8) | lsb[0]);
            }
            else
            {
                printf("Error reading file. %d bytes\n", read);
                exit(-1);
            }
        }
    }

    return 0;
}

struct HEADER readheader()
{
    int read = 0;
    unsigned char buffer4[4];
    unsigned char buffer2[2];

    struct HEADER header;
    // read header parts

    read = fread(header.riff, sizeof(header.riff), 1, stdin);
    if (DEBUG != 0)
        printf("(1-4): %s \n", header.riff);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    // convert little endian to big endian 4 byte int
    header.overall_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) |
                          (buffer4[3] << 24);
    if (DEBUG != 0)
        printf("(5-8) Overall size: bytes:%u, Kb:%u \n",
               header.overall_size,
               header.overall_size / 1024);

    read = fread(header.wave, sizeof(header.wave), 1, stdin);
    if (DEBUG != 0)
        printf("(9-12) Wave marker: %s\n", header.wave);

    read = fread(
        header.fmt_chunk_marker, sizeof(header.fmt_chunk_marker), 1, stdin);
    if (DEBUG != 0)
        printf("(13-16) Fmt marker: %s\n", header.fmt_chunk_marker);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    // convert little endian to big endian 4 byte integer
    header.length_of_fmt = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) |
                           (buffer4[3] << 24);
    if (DEBUG != 0)
        printf("(17-20) Length of Fmt header: %u \n", header.length_of_fmt);

    read = fread(buffer2, sizeof(buffer2), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u \n", buffer2[0], buffer2[1]);

    header.format_type = buffer2[0] | (buffer2[1] << 8);
    char format_name[10] = "";
    if (header.format_type == 1)
        strcpy(format_name, "PCM");
    else if (header.format_type == 6)
        strcpy(format_name, "A-law");
    else if (header.format_type == 7)
        strcpy(format_name, "Mu-law");
    if (DEBUG != 0)
        printf(
            "(21-22) Format type: %u %s \n", header.format_type, format_name);

    read = fread(buffer2, sizeof(buffer2), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u \n", buffer2[0], buffer2[1]);

    header.channels = buffer2[0] | (buffer2[1] << 8);
    if (DEBUG != 0)
        printf("(23-24) Channels: %u \n", header.channels);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    header.sample_rate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) |
                         (buffer4[3] << 24);

    if (DEBUG != 0)
        printf("(25-28) Sample rate: %u\n", header.sample_rate);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);

    header.byterate = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) |
                      (buffer4[3] << 24);
    if (DEBUG != 0)
        printf("(29-32) Byte Rate: %u , Bit Rate:%u\n",
               header.byterate,
               header.byterate * 8);
    read = fread(buffer2, sizeof(buffer2), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u \n", buffer2[0], buffer2[1]);

    header.block_align = buffer2[0] | (buffer2[1] << 8);
    if (DEBUG != 0)
        printf("(33-34) Block Alignment: %u \n", header.block_align);

    read = fread(buffer2, sizeof(buffer2), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u \n", buffer2[0], buffer2[1]);

    header.bits_per_sample = buffer2[0] | (buffer2[1] << 8);
    if (DEBUG != 0)
        printf("(35-36) Bits per sample: %u \n", header.bits_per_sample);

    read = fread(
        header.data_chunk_header, sizeof(header.data_chunk_header), 1, stdin);
    if (DEBUG != 0)
        printf("(37-40) Data Marker: %s \n", header.data_chunk_header);

    read = fread(buffer4, sizeof(buffer4), 1, stdin);
    if (DEBUG != 0)
        printf("%u %u %u %u\n", buffer4[0], buffer4[1], buffer4[2], buffer4[3]);
    if (read)
        read = 0;

    header.data_size = buffer4[0] | (buffer4[1] << 8) | (buffer4[2] << 16) |
                       (buffer4[3] << 24);
    if (DEBUG != 0)
        printf("(41-44) Size of data chunk: %u \n", header.data_size);

    // calculate no.of samples
    long num_samples =
        (8 * header.data_size) / (header.channels * header.bits_per_sample);
    if (DEBUG != 0)
        printf("Number of samples:%lu \n", num_samples);

    long size_of_each_sample = (header.channels * header.bits_per_sample) / 8;
    if (DEBUG != 0)
        printf("Size of each sample:%ld bytes\n", size_of_each_sample);

    // calculate duration of file
    float duration_in_seconds = (float)header.overall_size / header.byterate;
    if (DEBUG != 0)
        printf("Approx.Duration in seconds=%f\n", duration_in_seconds);

    return header;
}
