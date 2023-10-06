all: teeth capture test
teeth: teeth.c
	gcc teeth.c -o teeth -lfftw3 -lm

test: testfft.c libmylib.a 
	gcc -o testfft testfft.c -lasound -lm -lfftw3   -Wall -L. -lmylib

capture: capture.c defaultpulse.h libmylib.a 
	gcc -o capture capture.c -lasound -lm -lfftw3   -Wall -L. -lmylib
    
libmylib.o: mylib.c mylib.h

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o

clean:
	rm mylib.a mylib.o capture teeth libmylib.[oa] testfft
