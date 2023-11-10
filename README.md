# tico
mechanical watch analysis  on a Wav file
requires gnuplot and fftw3

```
$ sudo apt-get install fftw3-dev gnuplot
```

Example:

```
# make the wav file, use 48000Hz please
a=122;rm -f out.wav; (ffmpeg  -i Recorder/Voice\ $a.m4a out.wav) 2>/dev/null ; 
#help
tico -h
# use -t to get only the top of the sound
-r cuts off the last 10 seconds
-l cuts off the first 20 seconds
-d narrow search round max
-c crosscorrelation threshold
-v prints the gnuplot statement
-s split tick tock pulseshape
-e applies a gaussian window of s=8 

# generate a plot in lussen/www/tico.png and raw data
~/tico/tico/tico -l 4 -r 5 -s -e4 -d 50 < out.wav 

teeth -e4 -c 0.8 -l 4 -r 4 < out.wav
```
