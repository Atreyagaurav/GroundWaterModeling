set datafile separator ","

calibdata="./data/calib-trials-3.csv"

set yrange [0:100]
set multiplot layout 2,2
unset key
set grid ls 7 lc "gray"

set ylabel "MSE"
set xlabel "Kh"
plot calibdata using 1:4 with points ls 7 ps .4
set xlabel "RivKh"
plot calibdata using (($2+rand(0)/25)):4 with points ls 7 ps .4
set xlabel "Rech"
plot calibdata using (($3+rand(0)/5)):4 with points ls 7 ps .4
set xlabel "Trial #"
set yrange [-2:1]
set ylabel "NSE"
plot 1 with filledcurves y=0 fs transparent solid 0.4 fc "green",\
     .8 with lines ls 7 lc "blue",\
     calibdata using 0:5 with points ls 7 ps .4
