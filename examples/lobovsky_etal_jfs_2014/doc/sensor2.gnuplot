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
set key
set grid
# Set x parameters
set xlabel "t / (H/g)^0.5"
set xtic
set xrange[0:7]
# Set left y parameters
set ylabel "p / (rho g H)"
set ytic
set yrange[-0.2:2.5]
# Line styles
set style line 1 lt -1 lw 3
set style line 2 lt -1 lw 1 linecolor rgb "red"
set style line 3 lt -1 lw 1 linecolor rgb "blue"

# Plot
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$9/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'SPH' axes x1y1 with lines ls 1, \
     "doc/Fig30_filtered_1.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "doc/Fig30_filtered_2.dat" using 1:2 notitle axes x1y1 with lines ls 2, \
     "doc/Fig30_filtered_3.dat" using 1:2 title "Experiment 1, Wemmenhove et al. 2010" axes x1y1 with lines ls 3, \
     "doc/Fig30_filtered_4.dat" using 1:2 title "Experiment 2, Wemmenhove et al. 2010" axes x1y1 with lines ls 4




pause 5
replot
reread
