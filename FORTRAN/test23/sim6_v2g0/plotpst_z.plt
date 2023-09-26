set terminal postscript eps enhanced color font 'Helvetica,24'


set output "positions_z.eps" # Note, must write "setput" at first, otherwise there would be no data
set xlabel "T"
set ylabel "{/Symbol D}"
set title ""
set xrange [0:1e-1]
set xtics 2e-2
set mxtics 4
set yrange [0.0:3.0]
set ytics 0.5
set mytics 5

set label 1 sprintf("{/Symbol k}=%.1f",0.0) at 20,1.6 font "Arial, 24"
#set label 2 sprintf("{/Symbol D}=%.1f",1.0) at 15,1.6 font "Arial, 24"


dataf="k_0-d_0-positions.txt"

plot dataf using 1:2 index 1 title "" with lines lw 4 lc rgb 'black',\
#dataf using 1:2 index 2 title "" with lines lw 4 lc rgb '#a6cee3',\
#dataf using 1:2 index 3 title "" with lines lw 4 lc rgb '#1f78b4',\
#dataf using 1:2 index 4 title "" with lines lw 4 lc rgb '#b2df8a',\
#dataf using 1:2 index 5 title "" with lines lw 4 lc rgb '#33a02c',\
#dataf using 1:2 index 6 title "" with lines lw 4 lc rgb '#fb9a99',\
#dataf using 1:2 index 7 title "" with lines lw 4 lc rgb '#e31a1c',\
#dataf using 1:2 index 8 title "" with lines lw 4 lc rgb '#fdbf6f',\
#dataf using 1:2 index 9 title "" with lines lw 4 lc rgb '#ff7f00',\
#dataf using 1:2 index 10 title "" with lines lw 4 lc rgb '#cab2d6',\
#dataf using 1:2 index 11 title "" with lines lw 4 lc rgb '#6a3d9a',\
#dataf using 1:2 index 12 title "" with lines lw 4 lc rgb '#ffff99',\
#dataf using 1:2 index 13 title "" with lines lw 4 lc rgb '#b15928',\



set output
unset xrange
unset yrange
unset xtics
unset ytics
unset label