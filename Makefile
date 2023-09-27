teeth:
	gcc teeth.c -o teeth -lfftw3 -lm

capture:
	gcc -o capture capture.c -lasound -lm -lfftw3   -Wall 
