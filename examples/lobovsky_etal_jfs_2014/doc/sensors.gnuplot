#########################################################
#                                                       #
#    #    ##   #  #   #                           #     #
#   # #  #  #  #  #  # #                          #     #
#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #
#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #
#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #
#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #
#                            # #             #          #
#                          ##  #             #          #
#                                                       #
#########################################################
#
# plotSensors.gnuplot
#
# GnuPlot script that reads sensors file and plot in real time
# The pressure evolution at sensors.
#
#
#########################################################

# Set plot options
set multiplot
set key
set grid
set xlabel "t / (H/g)^0.5"
set xtic
set ylabel "p / (rho g H)"
set ytic
# Line styles
set style line 1 lt -1 lw 3
set style line 2 lt -1 lw 1 linecolor rgb "red"
set style line 3 lt -1 lw 1 linecolor rgb "blue"
set style line 4 lt -1 lw 1 linecolor rgb "magenta"

# Sensor 1
# ========
set xrange[0:7]
set yrange[-0.2:4.0]
set size nosquare 0.5,0.5
set origin 0.0,0.5
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$4/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'Sensor 1' axes x1y1 with lines ls 1, \
     "<awk '!/^($|[[:space:]]*#)/{t=$1/(0.3/9.81); p=100*$2; print t,p}' doc/Fig20_filtered_2.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "<awk '!/^($|[[:space:]]*#)/{t=$1/(0.3/9.81); p=100*$2; print t,p}' doc/Fig20_filtered_3.dat" using 1:2 notitle axes x1y1 with lines ls 2

# Sensor 2
# ========
set xrange[0:7]
set yrange[-0.2:2.5]
set size nosquare 0.5,0.5
set origin 0.5,0.5
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$9/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'Sensor 2' axes x1y1 with lines ls 1, \
     "doc/Fig30_filtered_1.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "doc/Fig30_filtered_2.dat" using 1:2 notitle axes x1y1 with lines ls 2, \
     "doc/Fig30_filtered_3.dat" using 1:2 title "Experiment 1, Wemmenhove et al. 2010" axes x1y1 with lines ls 3, \
     "doc/Fig30_filtered_4.dat" using 1:2 title "Experiment 2, Wemmenhove et al. 2010" axes x1y1 with lines ls 4

# Sensor 3
# ========
set xrange[0:7]
set yrange[-0.2:1.8]
set size nosquare 0.5,0.5
set origin 0.0,0.0
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$9/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'Sensor 3' axes x1y1 with lines ls 1, \
     "doc/Fig31_filtered_3.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "doc/Fig31_filtered_4.dat" using 1:2 notitle axes x1y1 with lines ls 2, \
     "doc/Fig31_filtered_2.dat" using 1:2 title "Experiment kleefsman et al. 2005" axes x1y1 with lines ls 3

# Sensor 4
# ========
set xrange[0:7]
set yrange[-0.2:1.2]
set size nosquare 0.5,0.5
set origin 0.5,0.0
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$19/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'Sensor 4' axes x1y1 with lines ls 1, \
     "doc/Fig29_filtered_3.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "doc/Fig29_filtered_4.dat" using 1:2 notitle axes x1y1 with lines ls 2, \
     "doc/Fig29_filtered_2.dat" using 1:2 title "Lee et al., 2002" axes x1y1 with lines ls 3

# Replot after 10 seconds
set nomultiplot
pause 10
replot
reread

