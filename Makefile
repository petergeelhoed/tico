FFTFLAGS= -lfftw3
SOUNDFLAGS= -lasound
MYLIBFLAGS= -L. -lmylib $(FFTFLAGS) $(SOUNDFLAGS) -pthread -Wpedantic -Wextra 
CFLAGS= -lm -Wall
CC= cc $(CFLAGS)

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
	$(CC) -o teeth teeth.c $(MYLIBFLAGS) 

testlinreg: testlinreg.c libmylib.a 
	$(CC) -o testlinreg testlinreg.c $(MYLIBFLAGS) 

fft: fft.c libmylib.a
	$(CC) -o fft fft.c $(MYLIBFLAGS) 

testfilter: testfilter.c libmylib.a 
	$(CC) -o testfilter testfilter.c  $(MYLIBFLAGS) 

testfft: testfft.c libmylib.a 
	$(CC) -o testfft testfft.c  $(MYLIBFLAGS) 

record: record.c libmylib.a 
	$(CC) -o record record.c $(MYLIBFLAGS) 

wav2raw: wav2raw.c 
	$(CC) -o  wav2raw wav2raw.c 
    
recali: recali.c libmylib.a 
	$(CC) -o recali recali.c $(MYLIBFLAGS) 
    
derivative: derivative.c libmylib.a 
	$(CC) -o derivative derivative.c  $(MYLIBFLAGS)
    
long: long.c defaultpulse.h libmylib.a 
	$(CC) -o long long.c $(MYLIBFLAGS) 
    
cap: cap.c libmylib.a 
	$(CC) -o cap cap.c $(MYLIBFLAGS)

capture: capture.c defaultpulse.h libmylib.a 
	$(CC) -o capture capture.c $(MYLIBFLAGS)

mylib.o: mylib.c mylib.h
	$(CC) -c  -o mylib.o mylib.c 

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o

tico:
	$(CC) -o tico tico.c  $(FFTFLAGS) $(SOUNDFLAGS) -pthread
