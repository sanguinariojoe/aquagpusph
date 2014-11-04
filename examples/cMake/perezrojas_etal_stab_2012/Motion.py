#########################################################################
#                                                                       #
#            #    ##   #  #   #                           #             #
#           # #  #  #  #  #  # #                          #             #
#          ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###           #
#          #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #          #
#          #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #          #
#          #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #          #
#                                    # #             #                  #
#                                  ##  #             #                  #
#                                                                       #
#########################################################################
#
#  This file is part of AQUA-gpusph, a free CFD program based on SPH.
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
#
#  AQUA-gpusph is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  AQUA-gpusph is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
#
#########################################################################

import math
import numpy as np
import aquagpusph as aqua


D = 0.062  # Tank width
g = 9.81   # Gravity acceleration [m/s2]
I0 = 26.9   # Polar moment of inertia [kg m2]
m = 4.978  # Moving mass weight [kg]
Sg = -29.2  # Satic moment of the rigid system [kg m]
Kdf = 0.540  # Dry friction damping coefficient 
Bphi = 0.326  # Linear damping coefficient

# Read the experimental data
data = []
f = open('@EXAMPLE_DEST_DIR@/Move/T_1-94_A100mm_water.dat', 'r')

# Discard header lines
line = f.readline()
line = f.readline()
n = 0
line = f.readline()
while(line):
    data.append([])
    words = line.strip().split(' ')
    for i in range(0,len(words)):
        data[n].append(float(words[i]))
    n = n+1
    line = f.readline()
f.close()
lineID = 0

# Starting mass motion data
Xi = data[0][1]
dXi = data[0][2]
# Starting Tank motion data
Theta = 0.0
dTheta = 0.0
ddTheta = 0.0
DT = 0.0

# Filter data
dThetaList = []
SMOOTH_WINDOW = 51

# Open the output file
f = open('Motion.dat', 'w')


def quit():
    """ Close the output file.
    """
    global f
    f.close()


def damping():
    """ Compute the tank structural damping moment.
    @return Tank structural damping moment.
    """
    global dTheta
    return -Kdf * np.sign(dTheta) - Bphi * dTheta


def angularForce(M):
    """ Compute the tank angular acceleration.
    @param M Measured fluid moment.
    @return Tank angular acceleration.
    """
    global Xi
    global dXi
    global Theta
    global dTheta
    Mdamp = damping()
    xi = Xi
    dxi = dXi
    K0 = I0 + m*xi*xi  # ddTheta term
    K1 = 2.0*m*xi*dxi  # dTheta term
    K2 = -g*Sg         # sin(Theta) term
    K3 = m*g*xi        # Independent term
    return (Mdamp + M - K1 * dTheta - K2 * math.sin(Theta) - K3 * math.cos(Theta)) / K0


def predictor():
    """ Performs predictor Leap-Frog stage.
    """
    global DT
    global Theta
    global dTheta
    global ddTheta
    dTheta = dTheta + DT * ddTheta;
    Theta = Theta  + DT * dTheta + 0.5 * DT * DT * ddTheta;


def corrector(dt, old):
    """ Performs predictor Leap-Frog stage.
    @param dt Time step
    @param old Previous ddTheta value
    """
    global DT
    global dTheta
    global ddTheta
    DT = dt
    dTheta = dTheta + 0.5 * DT * (ddTheta - old);


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    xx = np.asarray(x)

    if xx.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if xx.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return xx

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[xx[window_len-1:0:-1], xx, xx[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


def main():
    # Leap-frog Predictor
    predictor()
    # Get the new time and time step
    t = aqua.get("t")
    dt = aqua.get("dt")
    # Get new mass position
    global lineID
    global data
    global Xi
    global dXi
    while(t >= data[lineID][0]):
        lineID = lineID + 1
    factor = (data[lineID][0] - t) / (data[lineID][0] - data[lineID-1][0])
    Xi = (1.0 - factor)*data[lineID][1] + factor * data[lineID-1][1]
    dXi = (1.0 - factor)*data[lineID][2] + factor * data[lineID-1][2]
    # Get the experimentally obtained angle data
    exp_angle = (1.0 - factor) * data[lineID][4] + factor * data[lineID - 1][4]
    # Calculate ddTheta
    global ddTheta
    old = ddTheta
    M = 0.0  ###################################### FIXME
    ddTheta = angularForce(M)
    # Leap-frog Corrector
    corrector(dt, old)
    # Set the rotation angle. The angle is coming from a double integration
    # process, so it may be considered smooth enough.
    global Theta
    a = np.zeros(4, dtype=np.float32)
    a[0] = Theta
    aqua.set("motion_a", a)
    # Set the rotation velocity. The rotation velocity should be filtered due
    # to the large noise, which may ruin the simulation.
    global dTheta
    global dThetaList
    dThetaList.append(dTheta)
    dadt = np.zeros(4, dtype=np.float32)
    if len(dThetaList) < SMOOTH_WINDOW:
        dadt[0] = dTheta
    else:
        dadt[0] = smooth(dThetaList, SMOOTH_WINDOW)[-1]
        dThetaList[-SMOOTH_WINDOW:]
    aqua.set("motion_dadt", dadt)
    # Write output
    global f
    f.write('{}\t{}\t{}\t{}\t{}\n'.format(
        t, Xi, exp_angle,
        math.degrees(a[0]), math.degrees(dadt[0])))
    f.flush()
    return True
