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
	rm -f $(targets)aa


teeth: teeth.c libmylib.a libmysound.a
	$(CC) -o teeth teeth.c $(MYLIBFLAGS)  -lmysound

testlinreg: testlinreg.c libmylib.a  libmyfft.a
	$(CC) -o testlinreg testlinreg.c $(MYLIBFLAGS) -lmyfft

fft: fft.c libmylib.a libmyfft.a
	$(CC) -o fft fft.c $(MYLIBFLAGS) -lmyfft

testfilter: testfilter.c libmylib.a  libmyfft.a
	$(CC) -o testfilter testfilter.c  $(MYLIBFLAGS)  -lmyfft -lmylib

testfft: testfft.c libmylib.a libmyfft.a
	$(CC) -o testfft testfft.c  $(MYLIBFLAGS) -lmyfft -lmylib

record: record.c libmylib.a  libmysound.a libmyfft.a
	$(CC) -o record record.c $(MYLIBFLAGS)  -lmysound -lmyfft

wav2raw: wav2raw.c 
	$(CC) -o  wav2raw wav2raw.c 
    
recali: recali.c libmylib.a  libmysound.a libmyfft.a
	$(CC) -o recali recali.c $(MYLIBFLAGS)  -lmysound -lmyfft -lmylib
    
derivative: derivative.c libmylib.a 
	$(CC) -o derivative derivative.c  $(MYLIBFLAGS)
    
long: long.c defaultpulse.h libmylib.a libmysound.a libmyfft.a
	$(CC) -o long long.c $(MYLIBFLAGS)  -lmysound -lmyfft
    
cap: cap.c libmylib.a libmysound.a libmyfft.a
	$(CC) -o cap cap.c $(MYLIBFLAGS) -lmysound -lmyfft

capture: capture.c defaultpulse.h libmylib.a  libmysound.a libmyfft.a
	$(CC) -o capture capture.c $(MYLIBFLAGS) -lmysound -lmyfft

myfft.o: myfft.c myfft.h mylib.o
	$(CC) -c  -o myfft.o myfft.c -lmylib 

libmyfft.a: myfft.o myfft.h
	ar -rcs libmyfft.a myfft.o

mysound.o: mysound.c mysound.h
	$(CC) -c  -o mysound.o mysound.c 

libmysound.a: mysound.o mysound.h
	ar -rcs libmysound.a mysound.o

mylib.o: mylib.c mylib.h
	$(CC) -c  -o mylib.o mylib.c 

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o

tico:
	$(CC) -o tico tico.c  $(FFTFLAGS) $(SOUNDFLAGS) -pthread
