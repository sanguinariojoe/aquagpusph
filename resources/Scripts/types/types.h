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
 * @brief Generic types definition file.
 *
 * This file is just simply redirecting to either 2D.hcl or 3D.hcl, depending
 * on HAVE_3D
 */

#ifdef HAVE_3D
    #include "resources/Scripts/types/3D.h"
#else
    #include "resources/Scripts/types/2D.h"
#endif

#ifndef H2
#define H2 H*H
#endif

#ifndef SUPPORT2
#define SUPPORT2 SUPPORT*SUPPORT
#endif
