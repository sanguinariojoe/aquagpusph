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


D = 0.062  # Tank width
g = 9.81   # Gravity acceleration [m/s2]
I0 = 26.9   # Polar moment of inertia [kg m2]
m = 4.978  # Moving mass weight [kg]
Sg = -29.2  # Satic moment of the rigid system [kg m]
Kdf = 0.540  # Dry friction damping coefficient 
Bphi = 0.326  # Linear damping coefficient

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

# Starting mass movement data
Xi = data[0][1]
dXi = data[0][2]
# Starting Tank movement data
Theta = 0.0
dTheta = 0.0
ddTheta = 0.0
DT = 0.0
# Open output
try:
    f = open('Motion.dat', 'w')
except IOError:
    print('ERROR: Cannot open output data file.')
# File may still be opened when the simulation has finished

def init():
    """ Returns initial quaternion.
    @return [COR,X,Y], quaternion is defined by a center of rotation 
    and X,Y axes.
    """
    mCOR = [0.0,0.0,0.47]
    mX = [1.0,0.0,0.0]
    mY = [0.0,1.0,0.0]
    mZ = [0.0,0.0,1.0]
    return [mCOR, mX, mY, mZ]


def sign(value):
    """ Returns variable sign.
    @param value value to study
    @return value sign, if value is equal to 0, 0 will returned
    """
    if not value:
        return 0.0
    return math.copysign(1, value)


def damping():
    """ Returns tank damping moment.
    @return Tank damping moment.
    """
    global Kdf
    global dTheta
    return -Kdf * sign(dTheta) - Bphi * dTheta


def angularForce(Torque):
    """ Returns tank angular acceleration.
    @param Torque Measured fluid torque.
    @return Tank angular force.
    """
    global I0
    global m
    global g
    global Sg
    global Xi
    global dXi
    global Theta
    global dTheta
    Mdamp = damping()
    xi = Xi
    dxi = dXi
    K0 = I0 + m*xi*xi    # ddTheta term
    K1 = 2.0*m*xi*dxi    # dTheta term
    K2 = -g*Sg           # sin(Theta) term
    K3 = m*g*xi          # Independent term
    return (Mdamp - Torque - K1 * dTheta - K2 * math.sin(Theta) - K3 * math.cos(Theta)) / K0


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


def perform(COR, X, Y, Z, Torque, Force, t, dt):
    """ Returns quaternion evolution.
    @param COR Center of rotation at time t.
    @param X X axis at time t.
    @param Y Y axis at time t.
    @param Z Z axis at time t.
    @param Torque Measured fluid torque.
    @param Force Measured fluid force.
    @param t Time before integrate.
    @param dt Time step.
    @return [COR,X,Y], quaternion is defined by a center of rotation 
    and X,Y axes.
    """
    # Leap-frog Predictor
    predictor()
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
    # Get output angle data
    output = (1.0 - factor) * data[lineID][4] + factor * data[lineID - 1][4]
    # Calculate ddTheta
    global D
    global ddTheta
    old = ddTheta
    ddTheta = angularForce(Torque[0])
    # Leap-frog Corrector
    corrector(dt, old)
    # Write output
    global f
    global Theta
    f.write('{}\t{}\t{}\t{}\n'.format(t, Xi, math.degrees(Theta), output))
    f.flush()
    # Convert angle into quaternion
    mCOR = [0.0, 0.0, 0.47]
    mX   = [1.0, 0.0, 0.0]
    mY   = [0.0, math.cos(Theta), math.sin(Theta)]
    mZ   = [0.0, -math.sin(Theta), math.cos(Theta)]
    return [mCOR, mX, mY, mZ]
