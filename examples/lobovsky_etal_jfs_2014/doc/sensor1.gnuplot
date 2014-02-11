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
set yrange[-0.2:4.0]
# Line styles
set style line 1 lt -1 lw 3
set style line 2 lt -1 lw 1 linecolor rgb "red"

# Plot
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); p=$4/(997*9.81*0.3); print t,p}' Sensors.dat" using 1:2 title 'SPH' axes x1y1 with lines ls 1, \
     "<awk '!/^($|[[:space:]]*#)/{t=$1/(0.3/9.81); p=100*$2; print t,p}' doc/Fig20_filtered_2.dat" using 1:2 title "2.5% and 97.5% percentiles" axes x1y1 with lines ls 2, \
     "<awk '!/^($|[[:space:]]*#)/{t=$1/(0.3/9.81); p=100*$2; print t,p}' doc/Fig20_filtered_3.dat" using 1:2 notitle axes x1y1 with lines ls 2


pause 5
replot
reread
