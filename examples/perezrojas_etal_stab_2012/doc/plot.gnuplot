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

set grid
# Set x parameters
set xlabel "Time [s]"
set xtic
# Set left y parameters
set ylabel "Angle [deg]"
set ytic
# Line styles
set style line 1 lt -1 lw 1
set style line 2 lt -1 lw 1 lc rgb "#0000FF"
set style line 3 lt -1 lw 1 lc rgb "#FF0000"

# Plot
plot "zMovement.dat" using 1:3 title 'SPH' axes x1y1 with lines ls 1, \
     "zMovement.dat" using 1:4 title 'Experiments' axes x1y1 with lines ls 2
pause 5
replot
reread
