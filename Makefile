all: teeth capture
teeth: teeth.c
	gcc teeth.c -o teeth -lfftw3 -lm

capture: capture.c defaultpulse.h
	gcc -o capture capture.c -lasound -lm -lfftw3   -Wall 
