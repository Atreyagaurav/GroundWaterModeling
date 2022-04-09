set datafile separator ","

calibdata="../calib-trials-2.csv"

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
plot calibdata using 0:4 with points ls 7 ps .4

