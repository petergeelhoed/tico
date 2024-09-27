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


teeth: teeth.c libmylib.a libmysound.a
	$(CC) -o teeth teeth.c $(MYLIBFLAGS)  -lmysound

testlinreg: testlinreg.c libmylib.a 
	$(CC) -o testlinreg testlinreg.c $(MYLIBFLAGS) 

fft: fft.c libmylib.a
	$(CC) -o fft fft.c $(MYLIBFLAGS) 

testfilter: testfilter.c libmylib.a 
	$(CC) -o testfilter testfilter.c  $(MYLIBFLAGS) 

testfft: testfft.c libmylib.a 
	$(CC) -o testfft testfft.c  $(MYLIBFLAGS) 

record: record.c libmylib.a  libmysound.a
	$(CC) -o record record.c $(MYLIBFLAGS)  -lmysound

wav2raw: wav2raw.c 
	$(CC) -o  wav2raw wav2raw.c 
    
recali: recali.c libmylib.a  libmysound.a
	$(CC) -o recali recali.c $(MYLIBFLAGS)  -lmysound
    
derivative: derivative.c libmylib.a 
	$(CC) -o derivative derivative.c  $(MYLIBFLAGS)
    
long: long.c defaultpulse.h libmylib.a  libmysound.a
	$(CC) -o long long.c $(MYLIBFLAGS)  -lmysound
    
cap: cap.c libmylib.a libmysound.a
	$(CC) -o cap cap.c $(MYLIBFLAGS) -lmysound

capture: capture.c defaultpulse.h libmylib.a  libmysound.a
	$(CC) -o capture capture.c $(MYLIBFLAGS) -lmysound

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
