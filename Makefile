BIN=../bin
OBJ=../obj
LIB=../lib
SOUNDFLAGS= -L$(LIB) -lasound -lmysound
MYLIBFLAGS= -L$(LIB) -lmylib 
FFTFLAGS= -L$(LIB) -lfftw3 -lmyfft
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

all: $(targets)

clean:
	rm -f $(LIB)/* $(OBJ)/* $(BIN)/*

$(OBJ)/%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB)/lib%.a: $(OBJ)/%.o
	ar -rcs $@ $<

%: %.c $(LIB)/libmylib.a  $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/$@ $<  $(FFTFLAGS) $(SOUNDFLAGS) $(MYLIBFLAGS) 

capture: capture.c defaultpulse.h $(LIB)/libmylib.a  $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/capture capture.c $(MYLIBFLAGS) $(SOUNDFLAGS) $(FFTFLAGS)


testlinreg: testlinreg.c $(LIB)/libmylib.a  $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/testlinreg testlinreg.c $(MYLIBFLAGS) $(FFTFLAGS)

fft: fft.c $(LIB)/libmylib.a $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/fft fft.c $(MYLIBFLAGS) $(FFTFLAGS)

testfilter: testfilter.c $(LIB)/libmylib.a  $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/testfilter testfilter.c  $(MYLIBFLAGS)  $(FFTFLAGS) 

record: record.c $(LIB)/libmylib.a  $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/record record.c $(MYLIBFLAGS)  $(SOUNDFLAGS) $(FFTFLAGS)

recali: recali.c $(LIB)/libmylib.a  $(LIB)/libmysound.a $(LIB)/libmyfft.a
	$(CC) -o $(BIN)/recali recali.c  $(FFTFLAGS) $(MYLIBFLAGS)  $(SOUNDFLAGS)
    
derivative: derivative.c $(LIB)/libmylib.a 
	$(CC) -o $(BIN)/derivative derivative.c  $(MYLIBFLAGS)
    
tico: tico.c 
	$(CC) -o $(BIN)/tico tico.c  -lfftw3 -lasound -pthread

teeth: teeth.c $(LIB)/libmylib.a $(LIB)/libmysound.a
	$(CC) -o $(BIN)/teeth teeth.c $(MYLIBFLAGS)  -lmysound -lfftw3

wav2raw: wav2raw.c 
	$(CC) -o $(BIN)/wav2raw wav2raw.c 
    
