
#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#include "resources/Scripts/cfd/species/species_auxiliary.hcl"

#define E_ch 7000.0f
#define K_ch 1000000.0f

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