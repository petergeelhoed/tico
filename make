gcc tico.c  -lm -o tico -lfftw3 -I /usr/local/fftw/include -L /usr/local/fftw/lib
exit
a=122;rm -f out.wav; (ffmpeg  -i /home/charon/Voice\ Recorder/Voice\ $a.m4a out.wav) 2>/dev/null ; time ~/tico/tico/tico  -v -q 4000 -t   -l20 -r 10 -m5450 -s 100  -e8  -p15  < out.wav > raw
cat raw | plot 'u ($1/48):((1+9*($1/48<-8))*$2):($3):(sqrt($5**2)) w ye pal ;set xrange [-12:4]'


