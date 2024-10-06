BIN= ../bin
OBJ= ../obj
LIB= ../lib
SOUNDFLAGS= -lasound -lmysound 
MYLIBFLAGS= -lmylib 
MYSYNCFLAGS= -lmysync 
FFTFLAGS= -lfftw3 -lmyfft 
CFLAGS= -lm -Wall -pthread -Wpedantic -Wextra 
SANI_ADDR= -fsanitize=address -fno-omit-frame-pointer 
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
    tico\
    wav2raw

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
	rm -f $(LIB)/* $(OBJ)/* $(targets)

$(OBJ)/%.o: %.c %.h 
	$(CC) $(CFLAGS) -c -o $@ $< 

$(LIB)/lib%.a: $(OBJ)/%.o
	ar -rcs $@ $<

%: %.c $(LIB)/libmylib.a $(LIB)/libmysound.a $(LIB)/libmyfft.a $(LIB)/libmysync.a
	$(CC) -L$(LIB) -o $@ $<  $(SOUNDFLAGS) $(MYLIBFLAGS) $(FFTFLAGS) $(MYSYNCFLAGS)

#exectables
