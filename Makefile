# Directories
BIN = bin
OBJ = obj
LIB = lib

DIRS=\
     $(BIN)\
     $(OBJ)\
     $(LIB)

$(shell mkdir -p $(DIRS))

# Flags
CFLAGS = -g -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror # -Wconversion
# Uncomment the following line for optimized builds
# CFLAGS = -O3 -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -Wconversion
SANI_ADDR = -fsanitize=address -fno-omit-frame-pointer -static-libasan
SANI_THREAD = -fsanitize=thread -fno-omit-frame-pointer -static-libtsan
#CC= cc $(CFLAGS) $(SANI_ADDR)
#CC= cc $(CFLAGS) $(SANI_THREAD) 
CC = cc $(CFLAGS)

# Libraries
SOUNDFLAGS = -lasound -lmysound
MYLIBFLAGS = -lmylib
MYSIGNALFLAGS = -lmysignal
MYSYNCFLAGS = -lmysync
MYMATHFLAGS = -lmymath
FFTFLAGS = -lmyfft -lfftw3 

# Targets
TARGETS = \
    1dcor \
	derivative \
	capture \
	fft \
	ifft \
	recali \
	record \
	linregw \
	wav2raw

TESTTARGETS = \
	testfft \
	testfilter \
	testremove50 \
	testlinreg \
	testinvert \
	testmulmat \
	testmatlinreg \
	testmatlinregquad \
	testcalculateall \
	testcrosscorint \
	testtranspone 

# Outputs
OUTPUTS = \
	Fbase \
	capture \
	crosscor \
	filteredinput \
	flip \
	input \
	total

LIBS = \
	$(LIB)/libmylib.a \
	$(LIB)/libmysound.a \
	$(LIB)/libmysync.a \
	$(LIB)/libmysignal.a \
	$(LIB)/libmyfft.a \
	$(LIB)/libmymath.a

all: $(LIBS) $(TARGETS)
test: $(TESTTARGETS)

libs: $(LIBS)

install: $(TARGETS)
	mv $(TARGETS) $(BIN)

clean:
	rm -r $(LIB) $(OBJ) $(BIN) $(TARGETS) $(OUTPUTS) $(TESTTARGETS)

$(OBJ)/%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB)/libmylib.a: $(OBJ)/mylib.o $(OBJ)/crosscorint.o myarr.h
	ar -rcs $@ $< $(OBJ)/crosscorint.o

$(LIB)/lib%.a: $(OBJ)/%.o myarr.h
	ar -rcs $@ $<

%: %.c $(LIBS)
	$(CC) -L$(LIB) -o $@ $< $(SOUNDFLAGS) $(MYLIBFLAGS) $(FFTFLAGS) $(MYSYNCFLAGS) $(MYSIGNALFLAGS) $(MYMATHFLAGS) -lm  
