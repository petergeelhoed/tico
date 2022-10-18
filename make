gcc tico.c  -lm -o tico -lfftw3 -I /usr/local/fftw/include -L /usr/local/fftw/lib
exit
a=122;rm -f out.wav; (ffmpeg  -i /home/charon/Voice\ Recorder/Voice\ $a.m4a out.wav) 2>/dev/null ; time ~/tico/tico/tico  -v -q 4000 -t   -l20 -r 10 -m5450 -s 100  -e8  -p15  < out.wav > raw
cat raw | plot 'u ($1/48):((1+9*($1/48<-8))*$2):($3):(sqrt($5**2)) w ye pal ;set xrange [-12:4]'


gcc teeth.c  -lm -o teeth -lfftw3 -I /usr/local/fftw/include -L /usr/local/fftw/lib
 a=111; rm out.wav; ffmpeg -i ~/Voice/Voice\ $a.m4a  out.wav  &&  ./teeth l 5   -e 4   -d 50 -v -r 5  < ~/tico/github/tico/out.wav 


  #sort -k3,3g -k1,1g  raw |  plot 'u ($1-4000)/48:($2+100*$3):3  w l lw 2  pal    ; set xrange [-14:3]; set cbtics ("raw" 2,"abs(d/dt)" 4 , "gauss smooth" 6, "average pulseshape" 8 , "crosscor" 10 ); uns key ; set xlabel "time (ms)"; set ylabel "";  set format y; set colorbox user origin 0.01,0.33 size 0.02,0.53;'
