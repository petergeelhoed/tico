BIN= ../bin
OBJ= ../obj
LIB= ../lib
SOUNDFLAGS= -lasound -lmysound 
MYLIBFLAGS= -lmylib 
FFTFLAGS= -lfftw3 -lmyfft 
CFLAGS= -lm -Wall -pthread -Wpedantic -Wextra 
SANI_ADDR= -fsanitize=address -fno-omit-frame-pointer 
CC= cc $(CFLAGS)

targets=\
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

all:  $(libs) $(targets)

install: $(targets)
	mv $(targets) $(BIN)

libs= \
      $(LIB)/libmylib.a\
      $(LIB)/libmysound.a\
      $(LIB)/libmyfft.a 

clean:
	rm -f $(LIB)/* $(OBJ)/* $(targets)

$(OBJ)/%.o: %.c %.h 
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB)/lib%.a: $(OBJ)/%.o
	ar -rcs $@ $<

capture: capture.c $(LIB)/libmylib.a $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -L$(LIB) -o capture capture.c $(SOUNDFLAGS) $(FFTFLAGS) $(MYLIBFLAGS)

%: %.c $(LIB)/libmylib.a $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -L$(LIB) -o $@ $< $(SOUNDFLAGS) $(FFTFLAGS) $(MYLIBFLAGS)

#circular?
testlinreg: testlinreg.c $(LIB)/libmylib.a  $(LIB)/libmyfft.a
	$(CC) -L$(LIB) -o testlinreg testlinreg.c $(MYLIBFLAGS) $(FFTFLAGS)

fft: fft.c $(LIB)/libmylib.a $(LIB)/libmyfft.a
	$(CC) -L$(LIB) -o fft fft.c $(MYLIBFLAGS) $(FFTFLAGS)

tico: tico.c 
	$(CC) -o tico tico.c  -lfftw3 -lasound -pthread 

teeth: teeth.c
	$(CC) -o teeth teeth.c -lfftw3

wav2raw: wav2raw.c 
	$(CC)  -o wav2raw wav2raw.c 
    
derivative: derivative.c 
	$(CC) -o derivative derivative.c
    
