set terminal pdf size 20, 12
set output "data.pdf"

set xlabel 'turn[N]'
set ylabel 'position[mm]'
set grid

set multiplot layout 2,2

set yrange [-41.25:41.25]
set xtics 0, 10, 140

set title "Beam Position of #13"
plot './moment13.dat' using 0:3 with line linewidth 2 lt rgb '#4169e1' title 'x', \
     '' using 0:4 with line linewidth 2 lt rgb '#ff7f50' title 'y'

set title "Beam Position of #15"
plot './moment15.dat' using 0:3 with line linewidth 2 lt rgb '#4169e1' title 'x', \
     '' using 0:4 with line linewidth 2 lt rgb '#ff7f50' title 'y'

set yrange [0:400]
set xtics 0, 0.05, 1

set xlabel 'freqency'
set ylabel 'amplitude'

set title "Tune of #13"
plot './tune13.dat' using 1:2 with line linewidth 2 lt rgb '#4169e1' title 'x', \
     '' using 1:3 with line linewidth 2 lt rgb '#ff7f50' title 'y'

set title "Tune of #15"
plot './tune15.dat' using 1:2 with line linewidth 2 lt rgb '#4169e1' title 'x', \
     '' using 1:3 with line linewidth 2 lt rgb '#ff7f50' title 'y'

unset multiplot