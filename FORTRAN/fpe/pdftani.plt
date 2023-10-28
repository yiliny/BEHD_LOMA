set term gif size 800, 600
set term gif animate delay 2
set output "PDF_t.gif"
set xlabel "∆"
set ylabel "P"
set title ""
set key
set xrange [-0.01:3.0]
set xtics 0.5
set mxtics 4
ymax = 5.0e7
set yrange [-0.1:ymax]
set ytics ymax/10
set mytics 5

#GBf0="k_0-d_0-GBdist.txt"
GBf0="GBhdist.txt"
dataf0="k_0-d_0-PDF_ztime.txt"
dataf1="k_1-d_0-PDF_ztime.txt"
dataf2="k_2-d_0-PDF_ztime.txt"
dataf3="k_3-d_0-PDF_ztime.txt"
dataf4="k_4-d_0-PDF_ztime.txt"
dataf5="k_5-d_0-PDF_ztime.txt"
dataf6="k_6-d_0-PDF_ztime.txt"
dataf7="k_7-d_0-PDF_ztime.txt"
dataf8="k_8-d_0-PDF_ztime.txt"


do for [i=1:200-1]{
	#i=jj*100
	set label 1 sprintf("Time %i",i) at 0.5,ymax*0.7 font "Arial, 15"
	plot GBf0 using 1:($4*1.0) title "Boltzmann" with line dt 2 lt 6 lw 6 lc rgb 'grey',\
	dataf0 using 1:2 index i title "κ = 0.0e-1" with lines lw 2 lc rgb 'black',\
	dataf1 using 1:2 index i title "κ = 1.0e-1" with lines lw 2 lc rgb '#1f78b4',\
	dataf2 using 1:2 index i title "κ = 2.0e-1" with lines lw 2 lc rgb '#b15928',\
	dataf3 using 1:2 index i title "κ = 3.0e-1" with lines lw 2 lc rgb '#a65628',\
	dataf4 using 1:2 index i title "κ = 4.0e-1" with lines lw 2 lc rgb '#ff7f00',\
	dataf5 using 1:2 index i title "κ = 5.0e-1" with lines lw 2 lc rgb '#984ea3',\
	dataf6 using 1:2 index i title "κ = 6.0e-1" with lines lw 2 lc rgb '#4daf4a',\
	dataf7 using 1:2 index i title "κ = 7.0e-1" with lines lw 2 lc rgb '#377eb8',\
	dataf8 using 1:2 index i title "κ = 8.0e-1" with lines lw 2 lc rgb '#e41a1c',\
	
}

set output
unset xrange
unset yrange
unset xtics
unset ytics
unset label
unset key