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
set xlabel "t [s]"
set xtic
set xrange[0:5]
# Set left y parameters
set ylabel "P [Pa]"
set ytic
set yrange[-1000:5000]
# Line styles
set style line 1 lt -1 lw 1
set style line 2 lt  1 lw 1

# Plot
plot "<awk '{t=$1; p=$2; print t,100.0*p}' lateral_water_1x.txt" using 1:2 title 'Experimental' axes x1y1 with lines ls 2, \
     "<awk '{t=$1; p=$4+$9+$14+$19+$24+$29+$34+$39+$44;p=p/9.0; print t,p}' Sensors.dat" using 1:2 title 'SPH' axes x1y1 with lines ls 1

pause 5
replot
reread
