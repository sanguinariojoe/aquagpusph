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
import jonswap
import aquagpusph as aqua


H = {{DEPTH}}
G = {{G}}
HS = {{Hs}}
TP = {{Tp}}
SCALE_FACTOR = {{SCALE}}
NFREQS = 100
RAMP = 10.0


# Compute the spectra on real scale
T_bounds = [2.5, 18.0]
f = np.linspace(1 / T_bounds[1], 1 / T_bounds[0], num=NFREQS)
df = f[1] - f[0]
S = jonswap.jonswap(HS, TP, 2 * np.pi * f)

# Compute the waves set on real scale
Hw = 2.0 * np.sqrt(2.0 * S * df)

# Trasform the waves to model scale
Hw *= SCALE_FACTOR
f /= np.sqrt(SCALE_FACTOR)

# Compute the horizontal displacement of the flap at the free surface
# (far field Biesel transfer function for hinged-type wave-maker)
omega = 2 * np.pi * f
k = omega**2 / G
kh = k * H
# BUG: 2.0 here???
Hw_S0 = 4.0 * np.sinh(kh) * (1.0 - np.cosh(kh) + kh * np.sinh(kh)) / \
    (kh * (np.sinh(kh) * np.cosh(kh) + kh))
# Convert that to angles
S0 = Hw / Hw_S0
A = np.arctan(0.5 * S0 / H)
lag = np.random.rand(len(A)) * 2 * np.pi


def smoothstep(t):
   x = min(t / RAMP, 1)
   return x * x * (3 - 2 * x), 6 * x * (1 - x) / RAMP


def main():
    t = aqua.get("t")
    a = np.zeros(4, dtype=np.float32)
    ramp, drampdt = smoothstep(t)
    a[2] = ramp * np.sum(A * np.sin(omega * t + lag))
    aqua.set("motion_a", a)
    dadt = np.zeros(4, dtype=np.float32)
    dadt[2] = a[2] * drampdt + ramp * np.sum(
        A * omega * np.cos(omega * t + lag))
    aqua.set("motion_dadt", dadt)

    return True
