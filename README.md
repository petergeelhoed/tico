# tico
mechanical watch analysis  on a Wav file
requires gnuplot and fftw3
$ sudo apt-get fftw3-dev gnuplot

Example.
# make the wav file, use 48000Hz please
a=122;rm -f out.wav; (ffmpeg  -i Recorder/Voice\ $a.m4a out.wav) 2>/dev/null ; 
#help
tico -h
# use -t to get only the top of the sound
-r cuts off the last 10 seconds
-l cuts off the first 20 seconds
-v prints the gnuplot statement
-q moves the data in the modulo 8000 tics 
-m sets the center of the search area
-s the maximum distance to the -m
-e applies a gaussian window o f8 
-p 15 sets the number of teeth in the escapement wheel and prints the averaged data tot the stdout

# generate a plot in lussen/www/tico.png and raw data
~/tico/tico/tico  -v -q 4000 -t   -l20 -r 10 -m5450 -s 100  -e8  -p15  < out.wav > raw
#plot the raw data with a bit of a scale up for the first sound
cat raw | plot 'u ($1/48):((1+9*($1/48<-8))*$2):($3):(sqrt($5**2)) w ye pal ;set xrange [-12:4]'


