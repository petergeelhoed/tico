BIN= ../bin
OBJ= ../obj
LIB= ../lib
SOUNDFLAGS= -lasound -lmysound 
MYLIBFLAGS= -lmylib 
MYSYNCFLAGS= -lmysync 
FFTFLAGS= -lfftw3 -lmyfft 
CFLAGS= -g -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -Wconversion 
#CFLAGS= -O3 -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -Wconversion 
SANI_ADDR= -fsanitize=address -fno-omit-frame-pointer -static-libasan
#CC= cc $(CFLAGS) $(SANI_ADDR)
CC= cc $(CFLAGS)

targets=\
    derivative\
    capture\
    fft\
    recali\
    record\
    teeth\
    testfft\
    testfilter\
    testlinreg\
    testcalculateall\
    testcrosscorint\
    tico\
    wav2raw

outputs=\
        Fbase\
        capture\
        crosscor\
        filteredinput\
        flip\
        input\
        total


libs= \
      $(LIB)/libmylib.a\
      $(LIB)/libmysound.a\
      $(LIB)/libmysync.a\
      $(LIB)/libmyfft.a 

all:  $(libs) $(targets)
libs: $(libs)

install: $(targets)
	mv $(targets) $(BIN)

clean:
	rm -f $(LIB)/* $(OBJ)/* $(targets) $(outputs)

$(OBJ)/%.o: %.c %.h 
	$(CC) $(CFLAGS) -c -o $@ $< 

$(LIB)/libmylib.a: $(OBJ)/mylib.o $(OBJ)/crosscorint.o
	ar -rcs $@ $< $(OBJ)/crosscorint.o

$(LIB)/lib%.a: $(OBJ)/%.o
	ar -rcs $@ $<

%: %.c $(LIB)/libmylib.a $(LIB)/libmysound.a $(LIB)/libmyfft.a $(LIB)/libmysync.a
	$(CC) -L$(LIB) -o $@ $<  $(SOUNDFLAGS) $(MYLIBFLAGS) $(FFTFLAGS) $(MYSYNCFLAGS)

#exectables no conversion warnings.
tico: tico.c $(LIB)/libmyfft.a
	cc -g -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror  -L../lib -o tico tico.c   -lfftw3 -lmyfft

wav2raw: wav2raw.c $(LIB)/libmyfft.a
	cc -g -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror  -L../lib -o wav2raw wav2raw.c   -lfftw3 -lmyfft

teeth: teeth.c $(LIB)/libmyfft.a
	cc -g -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror  -L../lib -o teeth teeth.c   -lfftw3 -lmyfft
