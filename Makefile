SOUNDFLAGS= -L. -lasound -lmysound
MYLIBFLAGS= -L. -lmylib 
FFTFLAGS= -L. -lfftw3 -lmyfft
CFLAGS= -lm -Wall -pthread -Wpedantic -Wextra 
SANI_ADDR= -fsanitize=address -fno-omit-frame-pointer 
CC= cc $(CFLAGS)

targets=\
    cap\
    capture\
    derivative\
    fft\
    recali\
    record\
    teeth\
    testfft\
    testfilter\
    testlinreg\
    tico\
    wav2raw

libs=\
     mylib.a\
     myfft.a\
     mysound.aa
objects=\
        mylib.o\
        myfft.o\
        mysound.o

all: $(targets)

clean:
	rm -f $(targets) $(libs) $(objects)


testlinreg: testlinreg.c libmylib.a  libmyfft.a
	$(CC) -o testlinreg testlinreg.c $(MYLIBFLAGS) $(FFTFLAGS)

fft: fft.c libmylib.a libmyfft.a
	$(CC) -o fft fft.c $(MYLIBFLAGS) $(FFTFLAGS)

testfilter: testfilter.c libmylib.a  libmyfft.a
	$(CC) -o testfilter testfilter.c  $(MYLIBFLAGS)  $(FFTFLAGS) 

testfft: testfft.c libmylib.a libmyfft.a
	$(CC) -o testfft testfft.c  $(MYLIBFLAGS) $(FFTFLAGS)

record: record.c libmylib.a  libmysound.a libmyfft.a
	$(CC) -o record record.c $(MYLIBFLAGS)  $(SOUNDFLAGS) $(FFTFLAGS)

recali: recali.c libmylib.a  libmysound.a libmyfft.a
	$(CC) -o recali recali.c  $(FFTFLAGS) $(MYLIBFLAGS)  $(SOUNDFLAGS)
    
derivative: derivative.c libmylib.a 
	$(CC) -o derivative derivative.c  $(MYLIBFLAGS)
    
cap: cap.c libmylib.a libmysound.a libmyfft.a
	$(CC) -o cap cap.c $(MYLIBFLAGS) $(SOUNDFLAGS) $(FFTFLAGS)

capture: capture.c defaultpulse.h libmylib.a  libmysound.a libmyfft.a
	$(CC) -o capture capture.c $(MYLIBFLAGS) $(SOUNDFLAGS) $(FFTFLAGS)

myfft.o: myfft.c myfft.h mylib.o
	$(CC) -c  -o myfft.o myfft.c -lmylib

libmyfft.a: myfft.o myfft.h mylib.o
	ar -rcs libmyfft.a myfft.o mylib.o

mysound.o: mysound.c mysound.h
	$(CC) -c  -o mysound.o mysound.c 

libmysound.a: mysound.o mysound.h
	ar -rcs libmysound.a mysound.o

mylib.o: mylib.c mylib.h
	$(CC) -c  -o mylib.o mylib.c 

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o

tico: tico.c 
	$(CC) -o tico tico.c  -lfftw3 -lasound -pthread

teeth: teeth.c libmylib.a libmysound.a
	$(CC) -o teeth teeth.c $(MYLIBFLAGS)  -lmysound -lfftw3

wav2raw: wav2raw.c 
	$(CC) -o  wav2raw wav2raw.c 
    
