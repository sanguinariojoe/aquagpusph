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

# Set output file
set terminal postscript enhanced color
set output '| ps2pdf - images/lateral_water_1x_deleffe/press.pdf'

# Set plot options
set key
set grid
# Set x parameters
set xlabel "t [s]"
set xtic
set xrange[0:7]
# Set left y parameters
set ylabel "P [Pa]"
set ytic
set yrange[-1000:5000]
# Set right y parameters
set y2label "{/Symbol f} [deg]"
set y2tic
set y2range[-5:5]
# Line styles
set style line 1 lt -1 lw 1
set style line 2 lt  1 lw 2

# Plot
plot "<awk '{t=$1; phi=$3; print t,phi}' ../../examples/LateralWater_1x_DeLeffe/doc/lateral_water_1x.txt" using 1:2 title 'Roll angle' axes x1y2 with lines ls 1, \
     "<awk '{t=$1; p=$2; print t,100.0*p}' ../../examples/LateralWater_1x_DeLeffe/doc/lateral_water_1x.txt" using 1:2 title 'Pressure' axes x1y1 with lines ls 2
