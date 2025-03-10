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
FNAME = 'lateral_water_1x.txt'
T,_,A,DADT,_,_ = np.loadtxt(FNAME, delimiter='\t', skiprows=1, unpack=True)

def main():
    # Get the time instant
    t = aqua.get("t")
    # Interpolate the data
    a = np.zeros(4, dtype=np.float32)
    a[2] = np.radians(np.interp(t, T, A))
    dadt = np.zeros(4, dtype=np.float32)
    dadt[2] = np.radians(np.interp(t, T, DADT))
    # Send it to AQUAgpusph
    aqua.set("motion_a", a)
    aqua.set("motion_dadt", dadt)

    return True
