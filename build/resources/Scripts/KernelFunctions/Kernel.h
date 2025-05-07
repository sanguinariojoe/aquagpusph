/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Generic/automatic kernel header file.
 *
 * This file will read the definitions HAVE_3D and KERNEL_NAME, setting up and
 * including the appropiate file. In this way the user may easily select the
 * kernel to become applied using the modules at basic/kernels presets folder
 */

// Macro for adding quotes
#define KERNEL_STRINGIFY(X) KERNEL_STRINGIFY_AUX(X)
#define KERNEL_STRINGIFY_AUX(X) #X

// Concatenate macro
#define KERNEL_CAT(X,Y) KERNEL_CAT_AUX(X,Y)
#define KERNEL_CAT_AUX(X,Y) X##Y

// Kernel suffix
#ifdef HAVE_3D
    #define KERNEL_SUFIX 3D
#else
    #define KERNEL_SUFIX 2D
#endif

// Include the file
#define KERNEL_NAME_SUFIX KERNEL_CAT(KERNEL_NAME,KERNEL_SUFIX)
#include KERNEL_STRINGIFY(resources/Scripts/KernelFunctions/KERNEL_NAME_SUFIX.hcl)