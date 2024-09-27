CFLAGS= -pthread -lasound -lm -lfftw3 -Wall -L. -lmylib
CC= cc -c $(CFLAGS)

all: \
    cap\
    capture\
    derivative\
    fft\
    libmylib.a\
    recali\
    record\
    teeth\
    testfft\
    testfilter\
    testlinreg\
    tico\
    wav2raw

binaries=\
    mylib.[ao] \
    cap\
    capture\
    derivative\
    fft\
    libmylib.a\
    recali\
    record\
    teeth\
    testfft\
    testfilter\
    testlinreg\
    tico\
    wav2raw

clean:
	rm $(binaries)

teeth: teeth.c
	gcc teeth.c -o teeth -lfftw3 -lm

testlinreg: testlinreg.c libmylib.a 
	$(CC) -o testlinreg testlinreg.c

fft: fft.c
	cc $(CFLAGS) -o fft fft.c -lmylib

testfilter: testfilter.c libmylib.a 
	$(CC) -o testfilter testfilter.c 

testfft: testfft.c libmylib.a 
	$(CC) -o testfft testfft.c 

record: record.c libmylib.a 
	gcc -o record record.c  $(CFLAGS)

wav2raw: wav2raw.c 
	gcc -o wav2raw wav2raw.c 
    
recali: recali.c libmylib.a 
	$(CC) -o recali recali.c 
    
derivative: derivative.c
	gcc -o derivative derivative.c -lm -Wall 
    
long: long.c defaultpulse.h libmylib.a 
	gcc -o long long.c $(CFLAGS)
    
cap: cap.c libmylib.a 
	$(CC) -o cap cap.c -Wpedantic -Wextra 
    
capture: capture.c defaultpulse.h libmylib.a 
	gcc -o capture capture.c $(CFLAGS) 
    
libmylib.a: mylib.c mylib.h
	$(CC) -c -o mylib.o mylib.c 
	ar -rcs libmylib.a mylib.o
tico:
	$(CC) tico.c  -lm -o tico -lfftw3 -I /usr/local/fftw/include -L /usr/local/fftw/liba
