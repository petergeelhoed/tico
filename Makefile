# Directories
BIN = ../bin
OBJ = ../obj
LIB = ../lib

# Flags
CFLAGS = -g -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror # -Wconversion
# Uncomment the following line for optimized builds
# CFLAGS = -O3 -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -Wconversion
SANI_ADDR = -fsanitize=address -fno-omit-frame-pointer -static-libasan
#CC= cc $(CFLAGS) $(SANI_ADDR)
CC = cc $(CFLAGS)

# Libraries
SOUNDFLAGS = -lasound -lmysound
MYLIBFLAGS = -lmylib
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
	testfft \
	testfilter \
	testlinreg \
	testinvert \
	testmulmat \
	linregw \
	testmatlinreg \
	testmatlinregquad \
	testcalculateall \
	testcrosscorint \
	testtranspone \
	tico \
	wav2raw

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
	$(LIB)/libmyfft.a \
	$(LIB)/libmymath.a

all: $(LIBS) $(TARGETS)

libs: $(LIBS)

install: $(TARGETS)
	mv $(TARGETS) $(BIN)

clean:
	rm -f $(LIB)/* $(OBJ)/* $(TARGETS) $(OUTPUTS)

$(OBJ)/%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB)/libmylib.a: $(OBJ)/mylib.o $(OBJ)/crosscorint.o myarr.h
	ar -rcs $@ $< $(OBJ)/crosscorint.o

$(LIB)/lib%.a: $(OBJ)/%.o myarr.h
	ar -rcs $@ $<

%: %.c $(LIBS)
	$(CC) -L$(LIB) -o $@ $< $(SOUNDFLAGS) $(MYLIBFLAGS) $(FFTFLAGS) $(MYSYNCFLAGS) $(MYMATHFLAGS) -lm  

# Specific executables without conversion warnings
tico: tico.c $(LIB)/libmyfft.a
	cc -g -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -L$(LIB) -o tico tico.c -lfftw3 -lmyfft -lm 

wav2raw: wav2raw.c $(LIB)/libmyfft.a
	cc -g -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -L$(LIB) -o wav2raw wav2raw.c -lfftw3 -lmyfft -lm 

teeth: teeth.c $(LIB)/libmyfft.a
	cc -g -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -L$(LIB) -o teeth teeth.c -lfftw3 -lmyfft -lm 
