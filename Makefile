all: teeth capture testfft testlinreg testfilter tico record
teeth: teeth.c
	gcc teeth.c -o teeth -lfftw3 -lm

testlinreg: testlinreg.c libmylib.a 
	gcc -o testlinreg testlinreg.c -lasound -lm -lfftw3   -Wall -L. -lmylib

testfilter: testfilter.c libmylib.a 
	gcc -o testfilter testfilter.c -lasound -lm -lfftw3   -Wall -L. -lmylib

testfft: testfft.c libmylib.a 
	gcc -o testfft testfft.c -lasound -lm -lfftw3   -Wall -L. -lmylib

record: record.c libmylib.a 
	gcc -o record record.c -lasound -lm -lfftw3   -Wall -L. -lmylib
    
capture: capture.c defaultpulse.h libmylib.a 
	gcc -o capture capture.c -lasound -lm -lfftw3   -Wall -L. -lmylib
    
libmylib.o: mylib.c mylib.h

libmylib.a: mylib.o mylib.h
	ar -rcs libmylib.a mylib.o
tico:
	gcc tico.c  -lm -o tico -lfftw3 -I /usr/local/fftw/include -L /usr/local/fftw/liba
clean:
	rm  mylib.o capture teeth libmylib.[oa] testfft testlinreg testfilter
