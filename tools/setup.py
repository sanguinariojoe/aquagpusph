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

from distutils.core import setup

setup(
    name='AQUAgpusph-tools',
    version='1.4.6',
    author='J.L. Cercos Pita',
    author_email='jl.cercos@upm.es',
    packages=['aquagpusph_preprocessing',
              'aquagpusph_preprocessing.generator',
              'aquagpusph_preprocessing.mesh_loader'],
    scripts=['aquagpusph_preprocessing/AQUAgpusph-loadAbaqus',
             'aquagpusph_preprocessing/AQUAgpusph-loadGiD',
             'aquagpusph_postprocessing/pvd-locale'],
    url='http://canal.etsin.upm.es/aquagpusph',
    license='LICENSE',
    description='free SPH solver developed by the CEHINAV group.',
    long_description=open('../README.md').read(),
)

