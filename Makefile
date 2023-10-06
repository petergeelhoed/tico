all: teeth capture
teeth: teeth.c
	gcc teeth.c -o teeth -lfftw3 -lm

capture: capture.c defaultpulse.h libmylib.a 
	gcc -o capture capture.c -lasound -lm -lfftw3   -Wall -L. -lmylib
    
libmylib.o: mylib.c mylib.h

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o

clean:
	rm mylib.a mylib.o capture teeth libmylib.[oa]
