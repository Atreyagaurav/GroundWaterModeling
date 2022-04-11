set datafile separator ","
calibdata="./data/4_calibration_wells_model.csv"

set key top left

# from R
μ=-0.237
σ=1.4875

set yrange [0:]
set xrange [0:]

set linetype 1 lw 2 lc "blue" dt 1
set linetype 2 lw 1.5 lc "blue" dt 2
set linetype 3 lw 1 lc "blue" dt 3

plot calibdata using "h":"model_h" ls 7 title "Calib Data",\
     x title "1:1 Line" lt 1,\
     x-σ lt 2 notitle,\
     x+σ lt 2 title "± σ Line",\
     x-2*σ lt 3 notitle,\
     x+2*σ lt 3 title "± 2σ Line"

