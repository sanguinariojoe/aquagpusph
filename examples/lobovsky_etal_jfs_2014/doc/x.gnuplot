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
set key left top
# Set x parameters
set xlabel "t / (H/g)^0.5"
set xtic
set xrange[0:1.6]
# Set left y parameters
set ylabel "x / H"
set ytic
set yrange[0:1.6]
# Line styles
set style line 1 lt -1 lw 2
set style line 2 lt -1 lw 2 linecolor rgb "blue"

# Plot
plot "<awk '!/^($|[[:space:]]*#)/{t=$1/sqrt(0.3/9.81); x=-$2/0.3; print t,x}' Bounds.dat" using 1:2 title 'SPH' axes x1y1 with lines ls 1, \
     "doc/Fig12_01_ETSIN_Front_Dressler_Martin_2.dat" using 1:2 title "ETSIN experiment ... H=0.3m" axes x1y1 with lines ls 2

pause 5
replot
reread
