#********************************************************************************
#                                                                               *
#               *    **   *  *   *                           *                  *
#              * *  *  *  *  *  * *                          *                  *
#             ***** *  *  *  * *****  **  ***  *  *  ** ***  ***                *
#             *   * *  *  *  * *   * *  * *  * *  * *   *  * *  *               *
#             *   * *  *  *  * *   * *  * *  * *  *   * *  * *  *               *
#             *   *  ** *  **  *   *  *** ***   *** **  ***  *  *               *
#                                       * *             *                       *
#                                     **  *             *                       *
#                                                                               *
#********************************************************************************
#                                                                               *
#  This file is part of AQUAgpusph, a free CFD program based on SPH.            *
#  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>                 *
#                                                                               *
#  AQUAgpusph is free software: you can redistribute it and/or modify           *
#  it under the terms of the GNU General Public License as published by         *
#  the Free Software Foundation, either version 3 of the License, or            *
#  (at your option) any later version.                                          *
#                                                                               *
#  AQUAgpusph is distributed in the hope that it will be useful,                *
#  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
#  GNU General Public License for more details.                                 *
#                                                                               *
#  You should have received a copy of the GNU General Public License            *
#  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.          *
#                                                                               *
#********************************************************************************

from numpy import *

def length(v):
	""" Get vector length.
	@param v Vector
	@return length
	"""
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])