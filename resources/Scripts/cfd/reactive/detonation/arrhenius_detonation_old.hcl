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

#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#define E_ch 7000.0f
#define K_ch 1000000.0f
#define R_gas 8.31f

__constant float cps[4] = {14200.0f, 913.0f, 1040.0f, 2050.0f};
__constant float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
__constant float nus[4] = {-1.0f, -0.5f, 0.0f, 1.0f};
__constant float h_h20 = -285.83e3f/0.018f;
/*
float arrhenius_detonation(float zeta, float T);
float calc_cp_mix(float y_H2, float y_O2, float y_N2, float y_H2O);
void calc_gamma_cp_cv(float y_H2, float y_O2, float y_N2, float y_H2O, float* gamma, float* cp, float* cv);
float molar_mass_mixture(float y_H2, float y_O2, float y_N2, float y_H2O);
float M_zeta_dot(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O);
void w_rhos(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O, 
            float* w_rho_H2, float* w_rho_O2, float* w_rho_N2, float* w_rho_H2O, float* deintdt);
*/

float arrhenius_detonation(float zeta, float T) {

  return K_ch * (1.0f - zeta) * exp(-E_ch / T);
  
}

float molar_mass_mixture(float y_H2, float y_O2, float y_N2, float y_H2O){
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    return 1.0f /(y_H2 / Mis[0] + y_O2 / Mis[1] + y_N2 / Mis[2] + y_H2O / Mis[3]);

}

float calc_cp_mix(float y_H2, float y_O2, float y_N2, float y_H2O){

    //const float cps[4] = {14200.0f, 913.0f, 1040.0f, 2050.0f};
    return y_H2 * cps[0] + y_O2 * cps[1] + y_N2 * cps[2] + y_H2O * cps[3];

}

float M_zeta_dot(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O){
    float MMix;
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    
    float zx=MMix/Mis[0]*z;
    
    float zeta;
    if(zx>4.e-2f){    
        zeta = (z - y_H2) / (z+1.0e-12f);        
        return 1.0 / MMix * arrhenius_detonation(zeta, T);
    }else{
//        zeta = 1.0;
        return 0.0f;
    }

    

}

void X_from_Y(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* x_H2, __global float* x_O2, __global float* x_N2, __global float* x_H2O){
    float MMix;
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);

    *x_H2 = MMix / Mis[0] * y_H2;
    *x_O2 = MMix / Mis[1] * y_O2;
    *x_N2 = MMix / Mis[2] * y_N2;
    *x_H2O = MMix / Mis[3] * y_H2O;
    
    return;
}

void calc_gamma_cv(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* gamma, __global float* cv){

    float MMix, R_mix, cp_local, cv_local;
    
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    cp_local = calc_cp_mix(y_H2, y_O2, y_N2, y_H2O);
    
    R_mix = R_gas / MMix;

    //cv_local = cp_local / (cp_local-R_mix);
    cv_local = cp_local - R_mix;

    *gamma = cp_local / cv_local;
    
    *cv = cv_local;
    //printf("%f, %f, %f\n", MMix, cp_local, cv_local);

    return;
}

void calc_gamma_cp_cv(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* gamma, __global float* cv, __global float* cp){

    float MMix, R_mix, cp_local, cv_local;
    
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    cp_local = calc_cp_mix(y_H2, y_O2, y_N2, y_H2O);
    
    R_mix = R_gas / MMix;

    //cv_local = cp_local / (cp_local-R_mix);
    cv_local = cp_local - R_mix;

    *gamma = cp_local / cv_local;
    
    *cv = cv_local;
    //printf("%f, %f, %f\n", MMix, cp_local, cv_local);

    *cp = cp_local;
    
    return;
}


void w_rhos(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O, 
            __global float* w_rho_H2, __global float* w_rho_O2, 
            __global float* w_rho_N2, __global float* w_rho_H2O, 
            __global float* deintdt, __global float* zeta_dot){

    float M_zeta_dot_val = M_zeta_dot(z, T, y_H2, y_O2, y_N2, y_H2O);
    
    float MMix;
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    *zeta_dot = M_zeta_dot_val * MMix;

 
    //const float nus[4] = {-1.0f, -0.5f, 0.0f, 1.0f};
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    //const float h_h20 = -285.83e3f/0.018f;

    *w_rho_H2 = nus[0] * Mis[0] * M_zeta_dot_val;
    *w_rho_O2 = nus[1] * Mis[1] * M_zeta_dot_val;
    *w_rho_N2 = nus[2] * Mis[2] * M_zeta_dot_val;
    *w_rho_H2O = nus[3] * Mis[3] * M_zeta_dot_val;
    *deintdt -= h_h20 * *w_rho_H2O;
    
    //printf("%g\n",*w_rho_H2O);
    
    return;
}


#endif    // _ARRHENIUS_DETONATION_H_INCLUDED_