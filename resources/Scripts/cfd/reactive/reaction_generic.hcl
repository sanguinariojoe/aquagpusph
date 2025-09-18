
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

#ifndef _REACTION_GENERIC_H_INCLUDED_
#define _REACTION_GENERIC_H_INCLUDED_

#ifndef SPECIES_HEADER
#error "working with species requires to load a backend module"
#endif
#include SPECIES_HEADER

inline float give_zeta(float z, float y_0){

    // zeta prograss variable of combustion
    return (z - y_0) / (z+1.0e-12f);
}

void w_from_zeta_dot(__global species_t* w_rhos,
    __global float* deintdt, 
    __global float* zeta_dot, 
    float MMix){

        float inv_M_times_zeta_dot = *zeta_dot / MMix;
        
        const species_t Mis = MIS;
        const species_t nus = NUS;

        *w_rhos += nus * Mis * inv_M_times_zeta_dot;
      
        const species_t hplus_mass = HPLUS_MASS; 

        *deintdt -= dot(hplus_mass, *w_rhos);
        
        return;
    }

#endif    // _REACTION_GENERIC_H_INCLUDED_