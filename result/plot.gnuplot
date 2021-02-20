set terminal pdf size 20, 12
set output "data.pdf"

set xlabel 'Turn[N]'
set ylabel 'Sigma[mm*mm]'
set grid

set multiplot layout 2,2

set title "Beam Size"
plot './beamsize.dat' using 0:2 with points pointtype 2 lt rgb '#4169e1' title 'horizontal #13'
plot './beamsize.dat' using 0:3 with points pointtype 2 lt rgb '#ff7f50' title 'verticle #13'
plot './beamsize.dat' using 0:4 with points pointtype 2 lt rgb '#4169e1' title 'horizontal #15'
plot './beamsize.dat' using 0:5 with points pointtype 2 lt rgb '#ff7f50' title 'verticle #15'

unset multiplot