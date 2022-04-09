set datafile separator ","
calibdata="../calib-trials-3.csv"

set border 0
unset key
set xrange [] noextend
unset ytics
set style data parallelaxes

set palette model RGB functions gray*.8+.2,.7*(1-gray),0

# Turn on axis tics for the parallel axes
set for [i=1:5] paxis i tics

# Use the range commands to create an "interesting" plot.
# For suitable input data this is the sort of result you
# might get without axis range manipulation.

set paxis 1 range [*:120]
set paxis 2 range  [-30:*] reverse 
# set paxis 4 range  [-10:50]
# set paxis 5 range  [50:*] reverse
set paxis 5 tics left offset 4

plot calibdata using 4:4 lc palette title "RMSE",\
     '' using 5 title "NSE",\
     '' using 1 title "Kh",\
     '' using 2 title "RivKh",\
     '' using 3 title "Rech"

