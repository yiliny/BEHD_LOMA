set terminal postscript eps enhanced color font 'Helvetica,24'
set output "PDF_z.eps" # Note, must write "setput" at first, otherwise there would be no data
set xlabel "{/Symbol D}"
set ylabel "P_{eq}"
set title ""
set key
set xrange [-0.01:3.0]
set xtics 0.5
set mxtics 4
ymax = 8.0e7
set yrange [-0.1:ymax]
set ytics ymax/10
set mytics 5


#dataf="k_0-d_0-PDFz.txt"
GBf0="GBhdist.txt"
dataf0="k_0-d_0-PDF_ztime.txt"
dataf1="k_1-d_0-PDF_ztime.txt"
dataf2="k_2-d_0-PDF_ztime.txt"
dataf3="k_3-d_0-PDF_ztime.txt"
dataf4="k_4-d_0-PDF_ztime.txt"

i=200-1
set label 1 sprintf("Time %i",i) at 0.7,ymax*0.7 font "Arial, 15"
plot GBf0 using 1:($4*1) title "Boltzmann" with line dt 1 lt 6 lw 6 lc rgb 'black',\
dataf0 using 1:2 index i title "κ = 0.0e-1" with lines lw 4 lc rgb 'red',\
#dataf1 using 1:2 index i title "κ = 1.0e-1" with lines lw 2 lc rgb '#1f78b4',\
#dataf2 using 1:2 index i title "κ = 2.0e-1" with lines lw 2 lc rgb '#b15928',\
#dataf3 using 1:2 index i title "κ = 3.0e-1" with lines lw 2 lc rgb 'grey',\
#dataf4 using 1:2 index i title "κ = 4.0e-1" with lines lw 2 lc rgb '#ff7f00',\


set output
unset xrange
unset yrange
unset xtics
unset ytics
unset label