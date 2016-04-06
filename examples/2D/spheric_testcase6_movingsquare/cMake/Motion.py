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

import numpy as np
import os.path as path
import aquagpusph as aqua

# Read the experimental data
f = open(path.join('@EXAMPLE_DEST_DIR@/doc', 'Motion_Body.dat'), 'r')
lines = f.readlines()
f.close()
T = []
X = []
U = []
DUDT = []
for l in lines[2:]:
    l = l.strip()
    while l.find('  ') != -1:
        l = l.replace('  ', ' ')

    if l == '':
        continue
    t, dudt, u, x = map(float, l.split(' '))

    T.append(t)
    X.append(x)
    U.append(u)
    DUDT.append(dudt)
del f, lines

T = np.asarray(T)
X = np.asarray(X)
U = np.asarray(U)
DUDT = np.asarray(DUDT)

F = open('Motion.dat', 'w')

def main():
    # Get the time instant
    t = aqua.get("t")
    # Interpolate the data
    r = np.zeros(2, dtype=np.float32)
    r[0] = np.interp(t, T, X)
    u = np.zeros(2, dtype=np.float32)
    u[0] = np.interp(t, T, U)
    dudt = np.zeros(2, dtype=np.float32)
    dudt[0] = np.interp(t, T, DUDT)
    # Send it to AQUAgpusph
    aqua.set("motion_r", r)
    aqua.set("motion_drdt", u)
    aqua.set("motion_ddrddt", dudt)
    # Write output
    F.write('{}\t{}\t{}\t{}\n'.format(t, r[0], u[0], dudt[0]))
    F.flush()    
    return True