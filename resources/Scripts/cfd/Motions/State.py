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
import aquagpusph as aqua

motion_r = None
motion_a = None

def main():
    global motion_r, motion_a

    # Get the affected set of particles (it is also an identifier on the motion)
    motion_iset = aqua.get("motion_iset")

    # Check if they were not already set, eventually copying them from the
    # current motion state.
    if motion_r is None or motion_a is None:
        n_sets = aqua.get("n_sets")
        motion_r = [None] * n_sets
        motion_a = [None] * n_sets
    if motion_r[motion_iset] is None:
        motion_r[motion_iset] = aqua.get("motion_r")
    if motion_a[motion_iset] is None:
        motion_a[motion_iset] = aqua.get("motion_a")

    # Set the backuped state
    aqua.set("motion_r_in", motion_r[motion_iset])
    aqua.set("motion_a_in", motion_a[motion_iset])

    # Backuo the new state variables
    motion_r[motion_iset] = np.copy(aqua.get("motion_r"))
    motion_a[motion_iset] = np.copy(aqua.get("motion_a"))

    return True