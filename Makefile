BIN= ../bin
OBJ= ../obj
LIB= ../lib
SOUNDFLAGS= -lasound -lmysound 
MYLIBFLAGS= -lmylib 
MYSYNCFLAGS= -lmysync 
FFTFLAGS= -lfftw3 -lmyfft 
CFLAGS= -g -lm -Wall -pthread -Wpedantic -Wextra -Wsign-compare -Werror -Wconversion
SANI_ADDR= -fsanitize=address -fno-omit-frame-pointer 
CC= cc $(CFLAGS)

targets=\
    capture

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

