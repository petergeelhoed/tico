FFTFLAGS= -lfftw3
SOUNDFLAGS= -lasound
CFLAGS= -lm -Wall 
MYLIBFLAGS= -L. -lmylib
CC= cc -c $(CFLAGS)

targets=\
    mylib.o \
    libmylib.a \
    cap\
    capture\
    derivative\
    fft\
    long\
    recali\
    record\
    teeth\
    testfft\
    testfilter\
    testlinreg\
    tico\
    wav2raw

all: $(targets)

clean:
	rm -f $(targets)


teeth: teeth.c libmylib.a
	$(CC) teeth.c -o teeth 

testlinreg: testlinreg.c libmylib.a 
	$(CC) -o testlinreg testlinreg.c

fft: fft.c libmylib.a
	$(CC) -o fft fft.c 

testfilter: testfilter.c libmylib.a 
	$(CC) -o testfilter testfilter.c 

testfft: testfft.c libmylib.a 
	$(CC) -o testfft testfft.c 

record: record.c libmylib.a 
	$(CC) -o record record.c

wav2raw: wav2raw.c 
	$(CC) -o  wav2raw wav2raw.c 
    
recali: recali.c libmylib.a 
	$(CC) -o recali recali.c 
    
derivative: derivative.c
	$(CC) -o derivative derivative.c 
    
long: long.c defaultpulse.h libmylib.a 
	$(CC) -o long long.c
    
cap: cap.c libmylib.a 
	$(CC) -o cap cap.c -Wpedantic -Wextra 
    
capture: capture.c defaultpulse.h libmylib.a 
	$(CC) -o capture capture.c
    
mylib.o: mylib.c mylib.h
	$(CC) -o mylib.o mylib.c

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o
tico:
	$(CC) -o tico tico.c 
